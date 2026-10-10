# Implementation plan: hierarchical (F-model) allele-frequency prior

Branch `hier-prior`, from `master` (3.5.0). The sex-biased model (`dev`) and the
experimental `collapse-m` work are out of scope; they will be merged later.

## Motivation: the pooling artefact

BA3 gives each population an independent Dirichlet(1) prior on allele
frequencies. Each population's count table then pays its own cost for learning
the frequencies at every locus. When populations are weakly differentiated and
each table holds few gene copies (small samples, heavy missing data), moving
individuals' copies into another population's table saves more of that cost
than the Fst mismatch loses. The posterior then puts individuals into migrant
classes from a single source population, with m piled up near the 1/3 bound.
Below some amount of data per table pooling wins; above it, it loses, because
the mismatch cost grows linearly in the counts and the learning cost only
logarithmically.

Observed with simulated data under BA3's own model (6 pops x 15 indiv x 5,000
SNPs, Fst 0.05, `congen2026/sim/simulate_ba3.py`):

- No migration, F = 0, 30% missing genotypes: `-c` moves every individual in 5
  of 6 populations to second-generation migrants from one population by about
  iteration 4,500; the log-probability rises by about 19,000. Those populations
  have no age-0/1 individuals left, so their F is drawn from the prior
  (0.5 +/- 0.3), and m reaches about 0.26 on routes with no migration.
- Same data with 10% missing, or 0% missing: `-c` recovers m, the migrant
  routes and F correctly.
- The standard sampler does not leave its all-native start in 2M iterations
  (m about 0.017 everywhere, even on true routes with m = 0.1), so it does not
  show the mode; it is not evidence against it.
- The same artefact appeared with the standard sampler at Fst 0.005.

## Model

| Level | Current | Proposed |
|---|---|---|
| Ancestral frequency, locus l | none | pi_l ~ Dirichlet(1) |
| Population frequency | p_kl ~ Dirichlet(1), independent | p_kl given pi_l ~ Dirichlet(c_k pi_l), c_k = (1 - theta_k) / theta_k |
| Drift | none | theta_k in (0, 1), one per population |

theta_k -> 1 recovers the current model; theta_k -> 0 makes all populations
share pi_l. theta is a between-population drift (Fst-like) parameter, distinct
from BA3's within-population inbreeding F.

Unchanged: ancestry and its prior, m and its prior, F and the IBD indicators,
the age-2 phase marginal, missing-data handling, the Savage-Dickey test on m.

## Decisions

1. **theta per population** (Falush et al. 2003 F-model), not one shared theta.
2. **Prior on theta_k:** Uniform(0, 1) by default; a Beta prior concentrated on
   small values is tested in phase 6(c) and adopted only if justified.
3. **pi_l is sampled** by MH. Integrating pi_l out on a grid (biallelic loci
   only) is the fallback if it mixes poorly.
4. **Behind a flag** (`--hier`). Without it, output is byte-identical to 3.5.0.
   Whether it becomes the default is decided after validation.

## Collapsed sampler

Given pi_l and theta_k, p_kl is still Dirichlet, so `-c` integrates it out
exactly as now, with the Dirichlet-multinomial predictive

    P(next copy at (k, l) is allele a) = (n_kla + c_k pi_la) / (N_kl + c_k)

in place of (n_kla + 1) / (N_kl + A_l). pi_l and theta_k are sampled. The
marginal they need is the product over populations of

    Gamma(c_k) / Gamma(c_k + N_kl) * prod_a Gamma(c_k pi_la + n_kla) / Gamma(c_k pi_la)

which is computed from the existing count tables.

## Phases

**0. Confirm the premise.** Check that both samplers match the exact posterior
(`tests/check_exact.py`), including `exact_missing`. If `-c` passes, the pooled
mode is the real posterior of the current model and the plan proceeds. If it
fails, the `-c` bug is fixed first and this plan waits.

**1. Exact reference and prototype (Python, before any C++).**
- Extend `tests/exact_posterior.py` with the hierarchical prior: pi_l and
  theta_k integrated by quadrature on top of the existing enumeration. This is
  the exact reference for the later `check_exact.py` runs.
- Write a small Python collapsed sampler with the pi_l and theta_k moves; check
  it against the exact reference.
- Run the prototype on the pooling cases (no migration with 30% missing; Fst
  0.01) and on cases with true migration. Go/no-go: the pooled mode must
  disappear while true migration is still recovered.

**2. State and initialisation (C++).**
- `--hier` flag; pi[l][a], theta[k], and cached alpha[k][l][a] = c_k pi_la with
  per-(k, l) sums.
- Start pi_l at pooled sample frequencies and theta_k at a moment estimate of
  Fst between population k and the pooled sample.

**3. Generalise the collapsed hot path.**
- Real-valued alpha in `logAddRatioT`, `logGenoPredT`, `computeAddLogProb`,
  `addIndividual`, `removeIndividual`, `resampleIBD`, `collapsedLogLik`. The
  integer log table stays on the alpha = 1 path, so `-c` without `--hier` keeps
  its speed.
- The homozygote memo assumes one alpha for all loci; make it per locus or
  bypass it under `--hier`.
- Extend the `-d` self-check to the hierarchical marginal and run it with
  F > 0 as well as F = 0.

**4. New moves.**
- pi_l: MH per locus (logit random walk for biallelic loci, Dirichlet proposal
  centred on the current value for 3+ alleles), O(K A) from the count tables.
  Refresh one or a few random loci per iteration.
- theta_k: MH on log c_k, O(L A) per population, every n iterations; included
  in autotune.

**5. Standard sampler.** Add the hierarchical prior to the allele-frequency
move, and the same pi_l and theta_k moves using the explicit frequencies, so
that both samplers target the same posterior.

**6. Validation.**
- (a) `check_exact.py` against the phase-1 reference, both samplers, all
  `data/exact_*` sets.
- (b) Recovery with `simulate_ba3.py`, extended to draw frequencies from the
  hierarchical model with a set theta: no pooling at the 1/3 bound with no
  migration, across 0-50% missing and Fst 0.005-0.1; m and migrant classes
  recovered; theta_k recovers the simulated Fst.
- (c) Prior sensitivity on theta_k.
- (d) Empirical re-runs: the Malagasy palm and *Lumnitzera* SNP sets.
- (e) `tests/regress_default.sh`: without `--hier`, byte-identical to 3.5.0.
- (f) `tests/bench.py`: time and ESS per second against current `-c`.

**7. Output and documentation.** theta_k (posterior mean, SD) in the output and
the trace; optional pi output; manual section explaining the artefact; design
note `doc/hierarchical_prior.md` (maths and code map).

**Later (not on this branch).** Port to the sex-biased model on `dev` (X-linked
loci get their own pi_l; check hemizygous copies), and reconcile with
`collapse-m` (Rao-Blackwellized ancestry, stepping-stone, `--sdpool`). Keep pi
and alpha generic over ploidy so the port is mostly bookkeeping.

## Risks

- theta_k and migrant ancestry both explain populations looking alike; their
  mixing is the main thing phase 1 must show is acceptable.
- `--hier` loses the alpha = 1 integer-table shortcut; expect `-c --hier` to be
  1.5-3x slower per iteration than `-c`.
- This changes BA3's model, so it is a 3.6 or 4.0 change.
