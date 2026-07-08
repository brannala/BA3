# Implementation plan: integrate out allele frequencies (collapsed sampler)

Extends `doc/collapsed_allele_freqs.md` (math + code map) with a concrete,
phased build. The approach is validated by `sim/prototype_collapsed.py`
(marginal identity, incremental update, and collapsed==explicit posteriors with
inbreeding). Each phase is independently testable and leaves BA3 working; the
collapsed path lives behind a flag until it is validated, then becomes default.

## Architecture decision (fixed)

- **Persistent state:** per-individual ancestry (`migrantPopln`, `migrantAge`,
  `migrantSex` — unchanged) **plus** a per-individual, per-locus **IBD indicator
  `z`** for same-population diploid homozygous genotypes.
- **Sufficient statistic:** integer **count tables** `n[p][l][a]` (gene copies
  assigned to population `p`, locus `l`, allele `a`) and `N[p][l] = sum_a n`.
- Allele frequencies are **integrated out** (Dirichlet-multinomial); never
  sampled or stored. `F` becomes a **conjugate Beta-Bernoulli Gibbs** update from
  the `z` tallies.
- Concentration `alpha` (Dirichlet), default `1.0` (matches BA3's current
  implicit uniform prior). Single global constant `ALLELE_PRIOR_ALPHA`.

Why explicit `z` (not marginalized): the Dirichlet-multinomial requires integer
counts, and a homozygote contributes 1 copy (IBD) or 2 (outbred). `z` records
which, so the counts stay exact and updates stay O(1). `z` is Gibbs-sampled each
sweep; heterozygotes and age-2 hybrids and hemizygous males force `z` absent.

## New data structures

| Symbol | Type / size | Notes |
|---|---|---|
| `gCount` | `long[noPopln][noLoci][noAlleles[l]]` | same shape as the old `alleleFreqs`, integer |
| `gCountN` | `long[noPopln][noLoci]` | cached totals, avoids re-summing |
| `gIBD` | bit-packed `[noIndiv][ceil(noLoci/8)]` | `z` for same-pop diploid homozygotes; ~`noIndiv*noLoci/8` bytes |
| `ALLELE_PRIOR_ALPHA` | `const double` | Dirichlet concentration (1.0) |

Memory: `gCount` replaces `alleleFreqs`/`logAlleleFreqs`/`avgAlleleFreqs` (net
lower — one integer array instead of three double arrays). `gIBD` is small for
realistic sizes (e.g. 1000 indiv x 20k loci ~= 2.5 MB); bit-pack for the extreme
`MAXINDIV*MAXLOCI` corner.

## Core routines to add (all O(1) per gene copy)

```
logAddRatio(p,l,a)   = log( (n[p][l][a] + alpha) / (N[p][l] + A*alpha) )   // before incrementing
addCopy(p,l,a)       : n[p][l][a]++ ; N[p][l]++
removeCopy(p,l,a)    : n[p][l][a]-- ; N[p][l]--     // caller uses logAddRatio on the post-remove state for the reverse
```

Built on those:

```
removeIndividual(i)              // subtract i's gene copies from its current pop's cells (using i's ancestry + z)
computeAddLogProb(i, ancestry)   // log-prob of adding i's copies under a candidate ancestry, given current (without-i) counts;
                                 //   marginalizes z per homozygous locus via  F*pred1 + (1-F)*pred2  (prototype geno_pred)
addIndividual(i, ancestry)       // commit i's copies (sampling z per homozygous locus as it goes)
```

`computeAddLogProb` is the collapsed replacement for `logLik` in the ancestry
move: the MH genotype term becomes `computeAddLogProb(i, proposed) -
computeAddLogProb(i, current)` on the without-i counts (no cached per-individual
`logL`). Gene copies are tallied into the population each ancestry assigns them
to, including the age-2 male `sigma` (source vs native X) and hemizygous male X
(one copy) vs female X (two) — so the sex model composes unchanged.

## Phased build

**Phase 0 — scaffolding.** Add `ALLELE_PRIOR_ALPHA` and a flag (`gArgs.collapse`,
CLI `--collapse`, default off) so collapsed and current code coexist. No
behavior change.

**Phase 1 — count tables + full marginal.** Add `gCount`/`gCountN`; `initCounts`
(tally from the initial all-native state) replacing the role of
`getEmpiricalAlleleFreqs`; a `collapsedLogLik()` that sums the Dirichlet-
multinomial over all `(p,l)`. *Gate:* on a tiny fixed dataset, `collapsedLogLik`
equals the prototype's closed form.

**Phase 2 — incremental engine.** Implement `logAddRatio`/`addCopy`/`removeCopy`
and `removeIndividual`/`computeAddLogProb`/`addIndividual`. *Gate:* remove then
re-add an individual restores counts and marginal exactly; `computeAddLogProb`
matches recomputing `collapsedLogLik` from scratch (to 1e-9).

**Phase 3 — collapsed ancestry move (F = 0 first).** Rewire the ancestry
drop/add block (`main.cpp ~1211`, `dtLogL`) to use
`removeIndividual`/`computeAddLogProb`/`addIndividual` instead of `logLik`;
`fillMigrantCounts` accept path also updates `gCount`. Keep inbreeding off.
*Gate:* collapsed vs current sampler on a simulated `F=0` set (`simulate_sexbias.py`)
agree on `m` and per-individual ancestry posteriors within MC error.

**Phase 4 — inbreeding via IBD + F Gibbs.** Add `gIBD`; a per-iteration `z`
sweep (remove copy, resample `z` from `F*pred1 : (1-F)*pred2`, add copy) for
same-pop diploid homozygotes; replace the `FStat` MH block (`~1541`) with
`F_p ~ Beta(a0 + sum z_p, b0 + #diploid_p - sum z_p)`. *Gate:* collapsed vs
current on a simulated `F>0` set agree on `m` and `F` posteriors.

**Phase 5 — analytic frequency output.** *(Done — split from the deletion.)* The
standard sampler still carries the sex-biased model that the collapsed path lacks
until Phase 6, so its frequency MCMC and `alleleFreqs`/`avgAlleleFreqs` arrays are
**kept** for now. Instead, when `--collapse` is on, the sampling-step accumulator
draws `f` from its exact Dirichlet full-conditional (posterior parameters
`n[p][l][a]+alpha`) rather than reading the stale `alleleFreqs`; the reported
mean/SD then match the standard sampler's sampled-frequency output. *Gate (passed):*
collapsed vs standard allele-frequency output agree on both mean (max |Δ| < 0.01)
and SD. Deleting the frequency MCMC + arrays is deferred to Phase 6 (below), once
collapse is feature-complete and becomes the sole path.

**Phase 6 — sex model + finalize.** *(Sex integration done.)* The sex-biased model
now composes with `--collapse`: the count machinery was already sex-aware (age-2
male-X source/native copy assignment via `migrantSex`, hemizygous male X = one
copy), the `gamma` MH and `phi/rho` Gibbs blocks run unchanged, and the age-2
male-X `sigma` Gibbs gained a collapsed branch (remove the X copy, form the
log-odds from the count-based frequencies, resample, re-add — no `logAlleleFreqs`).
The `sexBiasModel`-under-`--collapse` force-off is removed. *Gate (passed):* on a
sex-biased set (2 pop, 300 indiv, 60 auto + 40 X, phiBreed 0.85) collapsed vs
standard agree on `m`, `F`, `rho`, `phiMove`, `phiBreed`, `gamma`, and the
Savage-Dickey BFs (all |z| < 0.25), collapsed recovers truth (score_recovery
PASS), and runs ~4x faster (9 s vs 40 s).

**Missing data:** *(decided — marginalize by omission.)* Under `--collapse` the
init-time random imputation of missing genotypes is skipped, so missing genotypes
stay `-1` and the count tables drop those gene copies (a fully-missing locus is
skipped; a partially-missing diploid contributes its one observed copy). This is
the exact Dirichlet-multinomial marginalization. `initCountsNative`,
`computeAddLogProb`, `addIndividual`, `removeIndividual`, `inbreedingUpdate`, and
the collapsed `sigma` Gibbs all guard `allele < 0`, and BA3's missing-data MH
(`if(!NOMISSINGDATA)`) stays gated off. The VCF reader marks a genotype missing
(both alleles `-1`) whenever *either* allele is absent, so genotype-level "partial"
missingness does not arise from VCF input. **Note:** the two paths therefore differ
on missing data — the standard sampler imputes missing genotypes from the current
(often native) population's frequencies, biasing migration downward, whereas
collapse marginalizes and stays close to the full-data estimate.

*Missing-data validation (both paths retained; run with matched seeds):*

- **Base model, 40 loci, 15% missing.** Before the fix (frozen random imputation)
  collapsed vs standard migration diverged, max |z| = 0.60. After the fix they
  agree within MC error. Anchored to the full-data estimate (m01 ~= 0.099, both
  paths agree): standard-missing ~= 0.077 (imputation biases it *down*),
  collapse-missing ~= 0.088 (near-unbiased, stable across 4 seeds).
- **Sex model, 100 loci (60 auto + 40 X), 15% missing (incl. missing male-X
  hemizygous loci).** Standard vs collapse agree within MC error on every
  parameter: `m` (|z| <= 0.14), `rho` (0.00), `phiMove` (0.04), `phiBreed` (0.34),
  `gamma` (0.50); collapse recovers truth (score_recovery PASS) and runs ~5.7x
  faster (10 s vs 58 s). Missing male-X copies are handled with no crash and
  `phiBreed` (informed by the X data) stays well-behaved.
- **Interpretation.** The imputation bias scales with the fraction of each
  individual's signal that is missing, so it is pronounced at 40 loci but nearly
  vanishes at 100 loci (each individual retains enough observed data). General
  statement: the paths are *not* identical on missing data (exact marginalization
  vs data-augmentation imputation); collapse is the correct treatment and is never
  worse; after the freeze-imputation fix both agree within MC error.

**Remaining (user decision).** Whether to (a) delete the now-redundant frequency
MCMC (`if(!NOALLELEMCMC)` block) + `alleleFreqs`/`logAlleleFreqs`/
`updateLogAlleleFreqs` and (b) flip `--collapse` to the default. Recommendation:
**keep both paths** — the standard sampler is a valuable cross-check and fallback,
and the collapse gates depend on running them side by side. Leave `--collapse`
opt-in until it has more field mileage. (The `z`-sweep is already thinned to every
`noIndiv` iterations — Phase 4.)

## Validation harness

- Unit gates in Phases 1-2 (marginal identity, incremental exactness).
- Posterior-equivalence gates in Phases 3-5: run collapsed and current samplers
  on the SAME `simulate_sexbias.py` dataset + seed and compare posteriors of
  `m`, `F`, `phi*`, `gamma`, and ancestry — they must agree within MC error
  (the prototype already shows this for a 2-pop mixture with inbreeding).
- Regression: existing autosomal-only and native-text inputs unchanged.

## Cost / risk

- **Per-iteration cost:** ancestry moves stay O(loci) with O(1) per-copy
  arithmetic (no `lgamma`); the `z` sweep is O(noIndiv*noLoci) simple arithmetic,
  replacing the current allele-freq MH sweep of similar order. Net expected win
  from deleting the frequency parameter block and improved mixing (ancestry
  decoupled from frequencies).
- **Risk:** this is a core rewrite; mitigated by the flag (both paths coexist)
  and the posterior-equivalence gates before the frequency MCMC is deleted.
- **Rollback:** until Phase 5, `--collapse off` is the untouched current sampler.

## Measured performance

Wall-clock, standard vs `--collapse`, matched seed and iteration count:

| Dataset | Indiv | Loci (X) | Iters | Standard | Collapse | Speedup |
|---|---|---|---|---|---|---|
| Simulated sex-biased | 300 | 100 (40) | 400k | 40 s | 9 s | ~4.3x |
| Simulated sex-biased, 15% missing | 300 | 100 (40) | 500k | 58 s | 10 s | ~5.7x |
| Brown bear (de Jong 2023) | 36 | 5102 (2582) | 200k | 51 s | 41 s | ~1.3x |

**The speedup scales with the individuals-to-loci ratio, not with size.** The
collapsed sampler's win comes from removing the allele-frequency parameters (and
their MCMC and the ancestry/frequency coupling), which dominates when individuals
outnumber loci. When loci greatly outnumber individuals (the bear set: 36 indiv,
5102 loci), per-iteration cost is instead dominated by O(loci) passes that *both*
samplers pay: the collapsed ancestry move makes ~4 loci-passes (remove + evaluate
current + evaluate proposed + re-add) to the standard move's ~1, and the age-2
male-X `sigma` Gibbs sweeps all X loci in both modes. Meanwhile the standard
frequency MCMC it eliminates is already cheap here (one locus per iteration, at
O(indiv) cost). So collapse still wins, but only modestly (~1.3x) in the
loci >> indiv regime, versus ~4-6x when indiv >= loci.

More precisely, the per-iteration advantage is a sum of terms, not purely O(indiv):
the removed frequency MCMC saves O(N), the removed `F`-statistic MCMC saves
O((N/P)*L) (it re-evaluates the chosen population's individuals across all loci),
and the collapsed ancestry move adds a ~3*O(L) penalty (extra loci-passes). Net
`~ a*N + b*(N/P)*L - c*L`: it grows with N and, via the F-stat term, with L as long
as N/P exceeds the ancestry penalty (~4). With `-u` (no inbreeding) the O((N/P)*L)
term drops out and the advantage reverts toward the pure-N frequency term.

### Scaling in N at fixed L

Why the speedup grows with N: the collapsed run time is *nearly* N-independent,
while the standard run time is ~linear in N. Base model (no sex), L = 100 loci,
2 populations, 150k iterations, varying total N:

| N | Collapse | Standard | std/col |
|---|---|---|---|
| 50 | 0.78 s | 2.16 s | 2.8x |
| 100 | 0.85 s | 5.07 s | 6.0x |
| 200 | 1.00 s | 11.9 s | 11.9x |
| 400 | 1.27 s | 26.6 s | 20.9x |

Collapse moves only 0.78 -> 1.27 s while N grows 8x: it fits `t ~= 0.71 + 0.0014*N`
-- a constant O(L) core (the one-individual ancestry move plus the inbreeding sweep
thinned to O(L) amortized) plus a *weak* linear-N tail from `fillMigrantCounts`
(O(N) per accepted ancestry move), sample-point sums over individuals, and (when
enabled) the O(N) phi/rho Gibbs tally. So the collapsed run time is **not strictly
independent of N, but only weakly dependent** (constant + small*N). The standard run
time is ~linear in N (its O(N) frequency MCMC and O((N/P)*L) F-stat MCMC), so the
ratio -- the speedup -- grows roughly linearly with N (2.8x at N=50 to 21x at
N=400). This is the same mechanism behind the dataset table above: the bear set
(N=36) sees a small speedup because standard's N-linear costs are still small there.

## Mixing (effective sample size)

Wall-clock speed is only half the story: the collapsed chain also *mixes* better,
because integrating out the allele frequencies removes the ancestry/migration <->
frequency coupling (a Rao-Blackwell / marginalization effect). Measured by the
integrated autocorrelation time (IAT) and effective sample size (ESS) of the
migration-rate trace columns (`sim/ess.py`, Geyer initial-positive-sequence
estimator; every-iteration trace, 80-100k post-burn-in samples, matched seed):

| Dataset | Standard median ESS/iter | Collapse median ESS/iter | Worst-column ESS (std -> col) |
|---|---|---|---|
| Simulated (300 indiv, 100 loci) | 0.0063 (IAT 167) | 0.0122 (IAT 82) | 398 -> 971 |
| Brown bear (36 indiv, 5102 loci) | 0.0043 (IAT 231) | 0.0067 (IAT 150) | **6 -> 184** |

Collapse mixes ~1.5-2x better per iteration in the median on both dataset shapes,
and the collapsed IATs are near-uniform across migration rates (e.g. 82.0-82.4 on
the simulated set) where the standard IATs vary (132-201). The sharpest effect is
on the *worst* coordinate: on the bear set the standard chain had a migration rate
that was essentially frozen (ESS 6 in 100k samples -- a small rate locked to
poorly-mixing frequencies), which collapse unsticks to ESS 184.

**Effective samples per second** = (ESS/iter) * (iters/sec) combines both factors
multiplicatively: ~1.9x mixing * ~4-6x speed ~= 8-11x on the simulated set, and
~1.5x mixing * ~1.3x speed ~= 2x on the bear set. So even in the loci >> indiv
regime where the raw speedup is modest, the mixing gain roughly doubles the real
efficiency, and the unstuck worst-parameter means far fewer iterations are needed
for a reliable estimate of every rate.

Reproduce: run each sampler with `-t -n 1` (dense trace), rename `BA3trace.txt`,
then `python sim/ess.py <label> <trace> <burnin>`.

## Open decisions to confirm before coding

1. `alpha` value (default 1.0 uniform; or a smaller symmetric value).
2. Missing-genotype handling under counts: keep BA3's imputation, or simply omit
   missing copies from the counts (cleaner — recommended).
3. `z`-sweep cadence: full sweep per iteration (simplest, correct) vs thinned.
4. Whether to keep both samplers long-term or make collapsed the only path.
