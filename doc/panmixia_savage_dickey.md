# Detecting over-split populations by model-selection on the number of demes

A test for **over-split populations** in BA3 input (near-identical demes labeled
as separate populations), which is the root cause of the gamma-identifiability
pathology (see the gamma/over-splitting notes): for near-identical populations
gamma is unidentifiable, migration inflates, and the sampler is seed-unstable.

## The right question, and the wrong ones

Compare **one panmictic population vs two populations _with migration_** — not vs
two _isolated_ populations. The isolated alternative is what several natural tests
implicitly assume, and it is wrong: under it, similar allele frequencies are an
improbable coincidence, so a Uniform frequency prior produces a Bartlett-Lindley
merge bias that grows with N and merges pairs up to F_ST ~ 0.5-0.9. Against the
_migration_ alternative, similar frequencies are expected (gene flow homogenizes),
so there is no similarity penalty.

Ruled out (verified), do not revisit:
- **Allele-count Dirichlet-multinomial BF** (shared f vs separate f): over-merges
  (Bartlett-Lindley; wrong, isolated-population alternative).
- **Genotype / Wahlund (F) term:** at within-group F ~ 0 the heterozygote factor
  `2^#het` is identical in the pooled and separate models, so it cancels in the BF
  and the genotype test reduces to the allele-count test. No extra power when F~0.
- **Savage-Dickey on an added Balding-Nichols divergence `d`:** the density ratio
  is cheap, but hosting it needs per-locus ancestral frequencies + `d` + a
  spike-and-slab -- nuisance parameters that also couple `f_A, f_B` and thereby
  break the collapsed sampler's independent-Dirichlet integration (its whole
  speed/mixing advantage). Not worth it.
- **Simulation-calibrated threshold:** the cutoff inherits the simulation's N,
  locus count and frequency spectrum -- an artifact, not a property of the data.

## The method: model-selection on the number of demes (stepping stone)

Do NOT pool a pair in isolation. Compare, on the **complete data**, the full
K-population migration model against the (K-1)-population model in which one pair is
**merged into a single population**:

```
log BF = log evidence(K-pop migration model)
       - log evidence((K-1)-pop model, one pair merged)
```

`BF > 1` keeps the finer split; `BF < 1` means over-split -> merge that pair. This
is model selection over population groupings: test the K-pop model against each
candidate merge of a suspect pair. It must be the joint model because merging A,B
changes their migration relationships with **every** other population (not just
A<->B), and it directly tests whether the A<->B migration structure is real or an
over-splitting artifact.

Both sides are full migration models, so **both need stepping stone** -- there is
no closed-form "one-population side" except the degenerate **K=2** case (merging
leaves a single panmictic population with no migration, whose evidence is the
Dirichlet-multinomial marginal integrated over one F).

### Why stepping stone is practical here

The collapsed sampler is the enabler, not incidental: with allele frequencies
integrated out analytically, each stepping-stone (power-posterior) run explores
only `(ancestry, migration rates, gamma, F)` -- a low-dimensional path that mixes
fast, so accurate thermodynamic integration needs fewer rungs and shorter runs even
on the full 13-population dace. Screen suspect pairs first (a quick F_ST scan, or a
short collapsed run flagging gamma-collapse or seed-instability) so only candidate
merges are tested, not all pairs.

## Plan

1. **Evidence estimator.** Add stepping-stone / power-posterior marginal-likelihood
   estimation on top of the collapsed MCMC (power the collapsed likelihood by beta
   from 0 to 1 across rungs). Validate the estimator on a case with a known/checkable
   evidence.
2. **K vs K-1 driver.** Given a candidate pair, build the (K-1)-population input
   (merge the two labels), run stepping stone for both models, report log BF.
3. **Screen + report.** Flag suspect pairs cheaply; run the test only on those;
   report per-candidate-merge log BF with a clear "populations X, Y are over-split /
   not separable" recommendation.

## Gate

On simulated data with known truth -- one population artificially split into two
labels vs two genuine demes (with and without migration, across a range of true
divergence) -- the log BF must merge the former and keep the latter. On the dace it
must flag AMA-DAMR / DWHS-LVD / DORB-RUP (F_ST 0.03-0.09 over-splits) while keeping
the F_ST 0.2+ pairs.

## Risk

For genuinely near-identical populations, "one population" and "two populations with
very high migration" can be **fundamentally confounded** (the same data). The joint
model + Occam (fewer demes) arbitrates, which is the correct behavior; the
simulation gate must confirm the evidence lands sensibly in that regime rather than
being prior-dominated. If it is prior-dominated there, "over-split / merge" is itself
the right answer for BA3, since recent-migration inference between such demes is
unreliable regardless.
