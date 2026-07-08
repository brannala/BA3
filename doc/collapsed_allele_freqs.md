# Collapsing (integrating out) population allele frequencies in BA3

**Status:** design scoping + validated prototype (`sim/prototype_collapsed.py`).

## Idea

BA3 currently **samples** allele frequencies `f_{p,l}` by MCMC. With BA3's
Dirichlet prior and a multinomial genotype likelihood (conditional on each gene
copy's ancestry), `f` can be **integrated out analytically** (Dirichlet-
multinomial / Polya). The sampler then works only in the ancestry states,
migration rates, sex parameters, and `F` — the entire `noPopln x noLoci x
noAlleles` frequency block disappears.

The prototype confirms this is exact (targets the same posterior), including the
inbreeding treatment:

```
[1] Dirichlet-multinomial marginal == MC integral            |diff| 0.005  OK
[2] incremental O(1) count update == recompute-from-scratch  exact         OK
[3] collapsed vs explicit Gibbs (2-pop mixture, inbreeding F)
      F   0.279 vs 0.277 ;  pi* 0.440 vs 0.439 ;  100% assignment agreement
```

## The collapsed likelihood

Maintain **count tables** `n[p][l][a]` = number of gene copies currently assigned
to population `p`, locus `l`, allele `a`, pooled over all individuals by their
ancestry (a native of `p`, an age-1 migrant from `p`, and the source half of an
age-2 hybrid from `p` all contribute). With a `Dirichlet(alpha,...,alpha)` prior
(`alpha = 1` reproduces BA3's current implicit uniform prior), the marginal
log-likelihood is a sum of Dirichlet-multinomial terms:

```
logL = sum_{p,l} [ lgamma(A*alpha) - lgamma(A*alpha + N_{p,l})
                   + sum_a ( lgamma(alpha + n_{p,l,a}) - lgamma(alpha) ) ]
```

where `A = noAlleles[l]`, `N_{p,l} = sum_a n_{p,l,a}`.

**Incremental updates are O(1)** (no `lgamma` calls). Adding one copy of allele
`a` to `(p,l)` multiplies the likelihood by `(n_{p,l,a} + alpha)/(N_{p,l} +
A*alpha)`; removing one divides by `(n-1+alpha)/(N-1+A*alpha)`. An ancestry move
that reassigns an individual moves its gene copies between two populations' count
cells at each locus — a handful of these ratios per locus. This replaces the
per-individual `logLik` product with a per-moved-copy ratio of the same cost.

## Inbreeding via IBD augmentation

The `F` term (`F*f_u + (1-F)*f_u*f_v` for a same-population homozygote) breaks the
plain multinomial. Fix it exactly by augmenting each **same-population diploid**
genotype with an IBD indicator `z ~ Bernoulli(F_p)`:

- `z = 1`: the two copies are IBD -> contributes **one** gene copy to `n`.
- `z = 0`: independent -> **two** copies.
- Heterozygotes force `z = 0`; age-2 hybrids are outbred (no `z`) — exactly BA3's
  current treatment (F applies to age 0/1, not age 2; never to hemizygous males).

Conditional on `z`, the counts are fixed and `f` integrates out as above.
`z_{i,l}` for homozygotes is Gibbs-sampled from its full conditional (prototype
`geno_pred` / the `pz1,pz0` block), and then **`F` is a conjugate Beta-Bernoulli
Gibbs update** from the `z` tallies:

```
F_p | . ~ Beta( a0 + #IBD-homozygotes_p , b0 + #non-IBD-diploid-genotypes_p )
```

replacing the current `F` Metropolis move.

## Mapping onto the current code

Current line numbers (branch `dev`), for reference:

| Current code | Change |
|---|---|
| `alleleFreqs` / `logAlleleFreqs` alloc (`~main.cpp:969-985`); `updateLogAlleleFreqs` (`~2962`) | **remove** — no sampled frequencies |
| `getEmpiricalAlleleFreqs` (`~2864`) | replace with **count-table init**: tally each individual's gene copies into `n[samplePopln][l][a]` under the initial all-native state |
| `logLik` / `oneLocusLogLik` (`~2591` / `~2745`) | replace with the count-based marginal; the ancestry move (`~1211`, `dtLogL`) uses the **incremental ratio** for the copies it moves instead of recomputing a per-individual product |
| allele-freq MCMC block `if(!NOALLELEMCMC)` (`~1440-1540`) | **delete**; add per-iteration Gibbs resampling of the IBD indicators `z` for homozygotes (updates `n` and the `F` tallies) |
| `FStat` MCMC block `if(!NOFSTATMCMC)` (`~1541`) | replace with the **Beta-Bernoulli Gibbs** update of `F_p` from `z` |
| freq averaging (`~1987-1990`) and freq output (`~2234-2260`) | report the **analytic** posterior of `f` on demand: posterior mean `(n_{p,l,a}+alpha)/(N_{p,l}+A*alpha)` (and the Dirichlet SD), accumulated over samples — no `avgAlleleFreqs` array needed |
| ancestry accept updates `fillMigrantCounts` etc. | additionally move the individual's gene copies between count cells (and keep the `F` tallies consistent) |

New global state: `long n[noPopln][noLoci][maxAlleles]` (counts) and `N[noPopln][noLoci]`; per-individual per-homozygous-locus `z` (can be packed, or recomputed since it only exists for homozygotes). Memory is comparable to the old frequency arrays but integer, and there is no longer a per-locus MH sweep over frequencies.

## Compatibility with the sex-biased model

Clean: a gene copy is tallied into whatever population its ancestry assigns it to,
including the age-2 male `sigma` (source vs native X). Hemizygous male X
contributes one copy, female X two — as already coded. `m`, `phiMove`,
`phiBreed`, `rho`, `gamma` do not touch `f`, so they are unaffected. Missing
genotypes: impute as now, or (better) simply omit missing copies from the counts.

## Benefits

- **Speed / memory:** eliminates the dominant per-iteration cost (the allele-freq
  MH sweep) and the frequency-array storage on large SNP panels; ancestry moves
  stay O(loci) with O(1) per-copy arithmetic.
- **Mixing:** collapsing decouples ancestry states from frequencies (Rao-Blackwell
  / collapsed-Gibbs) — directly improves the ancestry moves that are the crux.
- **Cleaner `F`:** conjugate Gibbs instead of MH.
- **Same posterior** (a marginalization, not a model change); allele-frequency
  posteriors remain reportable analytically.

## Validation plan

1. `sim/prototype_collapsed.py` — done (marginal identity, incremental update,
   collapsed==explicit with inbreeding).
2. In BA3: build the collapsed likelihood behind a flag, run collapsed vs current
   on a simulated dataset (`sim/simulate_sexbias.py`), and confirm the posteriors
   of `m`, `F`, `phi*`, `gamma`, and per-individual ancestry agree within MC
   error. Then delete the frequency MCMC.

## Effort / risk

Medium-large: a core rewrite of the likelihood + allele-freq/F machinery
(count tables, incremental ratios, IBD/`F` Gibbs, analytic frequency output).
Best done as its own phase on a branch, validated against the current sampler
before the frequency MCMC is removed. Confirm the target concentration `alpha`
(current implicit prior is uniform, `alpha = 1`).
