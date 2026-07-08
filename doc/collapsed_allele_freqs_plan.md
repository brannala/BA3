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

**Phase 5 — remove the frequency MCMC + analytic output.** Delete the
`if(!NOALLELEMCMC)` block (`~1440-1540`) and the `alleleFreqs`/`logAlleleFreqs`/
`avgAlleleFreqs` arrays and `updateLogAlleleFreqs`. Replace frequency averaging/
output with the analytic Dirichlet posterior mean `(n[p][l][a]+alpha)/(N[p][l]+A*alpha)`
(and its SD), accumulated at sampling steps. *Gate:* `score_recovery.py` still
passes on the `phi_move/phi_breed/gamma` simulations.

**Phase 6 — finalize.** X/sex + missing-data integration checks (impute or omit
missing copies from counts — decide and document); optional `z`-sweep thinning if
the O(noIndiv*noLoci) sweep dominates on large panels; flip `--collapse` to
default (or keep both paths).

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

## Open decisions to confirm before coding

1. `alpha` value (default 1.0 uniform; or a smaller symmetric value).
2. Missing-genotype handling under counts: keep BA3's imputation, or simply omit
   missing copies from the counts (cleaner — recommended).
3. `z`-sweep cadence: full sweep per iteration (simplest, correct) vs thinned.
4. Whether to keep both samplers long-term or make collapsed the only path.
