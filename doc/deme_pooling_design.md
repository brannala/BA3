# Testing whether sampled populations are one deme: design note

Status: design, 2026-10-10. Branch `hier-prior`. Nothing implemented.
Plan: section 8 (fixed pair, then any pair, then the full partition model).

## 1. Problem

BA3 treats every sampling location as a separate population. When two
locations, A and B, actually sample one panmictic population ("oversplitting"),
or differ only slightly, an individual sampled in B fits equally well as a
resident of B or as a migrant from A. Its ancestry, and so m_BA and m_AB, are
then determined by the prior rather than the data. The flat allele-frequency
prior amplifies this into the pooling artefact (whole populations assigned as
migrants, m near its bound). It would be useful to have a statistically
rigorous criterion for treating A and B as one population.

Approaches already tried or rejected:

- **Stepping-stone marginal likelihoods (`-S`, `collapse-m`)**: estimates the
  whole-data evidence of pooled and split models separately. On large SNP data,
  treating linked SNPs as independent (a composite likelihood) and multimodal
  ancestry gave a seed-to-seed SD of up to ~1,200 nats in log p(G), larger than
  the differences being tested.
- **Savage-Dickey on equal frequencies (`-P`, `collapse-m`)**: a posterior mean
  of a product of per-locus terms; its variance is unbounded with many loci.
- **Migrant detectability** (expected per-individual log-likelihood ratio,
  summed over loci): confounds population differentiation with the information
  in one individual's multilocus genotype. It grows with the number of loci at
  a fixed Fst, and is biased upward under panmixia by noise in the frequency
  estimates.
- **A population tree with nested Dirichlet frequencies** (shared node for A
  and B): a nested Dirichlet matches only the second moments of drift. It is
  not the Wright-Fisher process (no loss or fixation, wrong shape), the leaves'
  joint distribution is not that of a tree, and internal-node frequencies do
  not integrate out. Proper tree models (SNAPP-type coalescent likelihoods,
  diffusion, Gaussian approximations as in TreeMix) are different machinery.

## 2. Model: populations grouped into demes

The single-level hierarchical prior (`-H`) has a population-genetic
justification as Wright's island model. At migration-drift equilibrium a
deme's allele frequencies have, in the diffusion limit, a Dirichlet
distribution around the migrant-pool frequencies pi, with p ~ Dir(c pi),
c = 4Nm = (1-theta)/theta. The proposal keeps this model and changes only which
samples belong to which deme.

A **partition** g assigns each sampled population s to a deme g(s). Each deme
d has

- allele frequencies p_d ~ Dir(c_d pi) per locus (`-H`), or Dir(1) (flat
  prior), with theta_d ~ U(0, 1);
- an inbreeding coefficient F_d ~ U(0, 1);
- a migration row m_d over the demes, flat on {m_dj >= 0, sum_{j!=d} m_dj <= 1/3}.

An individual sampled in population s is a resident of deme g(s), or a first-
or second-generation migrant from another deme, with BA3's usual probabilities
1-3 sum_j m, m and 2m. The split model (every population its own deme) is
BA3 as it is now.

**Pair test.** M1 (split): A and B are separate demes. M0 (merged): A and B are
one deme with one frequency vector per locus, one F, one theta and one
migration row. All other demes are unchanged, except that their migration rows
now have one source, A∪B, in place of two.

**General form.** Posterior over partitions of the K sampled populations,
reported as posterior probabilities of partitions and pairwise co-membership
probabilities. Demes are labelled by their sets of sampled populations, so
there is no label switching.

Both M0 and M1 are proper island models. The test asks whether A and B have
identical allele frequencies (exact panmixia between them), not whether they
form a clade.

**Relation to `-H`.** `-H` is the special case in which every sampled
population is its own deme: one draw p_q ~ Dir(c_q pi) and one theta_q per
sampled population, with pi shared. The extension changes only the map from
samples to draws. A deme's frequencies are drawn once around the shared pi,
and every sampled population in the deme takes its individuals from that same
draw. Wrongly splitting a deme means taking two multinomial samples from one
frequency vector, with no between-sample variance: the Dirichlet-multinomial
with theta = 0 between the two samples. The star prior of `-H` cannot express
that, because it draws the two samples' frequencies independently around pi,
so the split model can make them identical only through migration. In
simulations, with 4 demes, Fst 0.05, 5,000 SNPs and one deme's 100 individuals
split at random into two labels of 50, `-c` and `-c -H` both put m between the
halves at the bound (0.30-0.31) in one direction and 0.21-0.31 in the other,
reported decisive Bayes factors for migration in both directions, and often
returned F at its prior (0.5, SD 0.29) for one half.

A continuous alternative, a within-deme drift theta_w with p_A, p_B ~
Dir(c_w p_AB) and pooling at theta_w = 0, is not used. It needs a second
Dirichlet level, which only approximates drift (see section 1), while the
partition model keeps a single draw around pi for every deme.

## 3. Why a merge/split move is tractable

Under the collapsed sampler, almost everything that changes dimension between
M0 and M1 integrates out in closed form, given the latent ancestry and IBD
indicators:

| Parameter | Integrated how | Effect of merging A and B |
|---|---|---|
| Allele frequencies p | Dirichlet-multinomial (`-c`), with c_d pi under `-H` | add the A and B count tables |
| F | Beta-binomial on IBD indicators: p(Z \| T) = B(1+Z, 1+T-Z) | add the IBD totals Z and T |
| m | Polya urn on the ancestry counts (u = 3m is flat Dirichlet) | merge the A and B categories in every row; merge the A and B rows |
| theta | not closed form: 1-D quadrature or proposal (below) | one theta in place of two |
| pi | shared by all demes | unchanged |

The posterior probability of M0 is then estimated from one chain that moves
between M0 and M1. No whole-data marginal likelihood is needed, which avoids
the stepping-stone failure.

**m handling.** Integrate m out only within the move: accept on the
m-marginal, then immediately redraw m from its full conditional (a Dirichlet
given the new ancestry counts). This is a valid partially collapsed Gibbs step
(van Dyk & Park 2008), provided the m draw immediately follows the move. It
keeps the rest of the `-c` sampler unchanged and does not require `-C`.

**theta handling.** Either (a) integrate theta_AB numerically in the
acceptance ratio (1-D Gauss-Legendre, O(L x nodes) per proposal), or (b)
propose the new thetas from an independence proposal fitted from the counts
(a moment estimate) and include the proposal densities in the Hastings ratio.
Option (a) is simpler to get right; (b) is cheaper. To be decided after
profiling.

## 4. Ancestry mapping

**Merge (M1 -> M0)** is deterministic except for IBD indicators.

| Individual, sampled in | Ancestry under M1 | Under M0 |
|---|---|---|
| A (or B) | resident | resident of A∪B |
| A | 1st-gen from B | resident of A∪B (both copies from A∪B; IBD indicator kept, counted toward F_AB) |
| A | 2nd-gen from B | resident of A∪B: **needs IBD indicators** at homozygous loci (age 2 has none in BA3); propose from Bernoulli(F_AB) at the current F_AB, or from their conditional |
| A | from C (any age) | from C into A∪B |
| C | from A or B (any age) | from A∪B into C |

**Split (M0 -> M1)** must propose:

- for each individual sampled in A or B and resident in A∪B: resident, or a
  1st- or 2nd-gen migrant from the other member of the pair. The proposal is
  sequential, drawing each individual from its conditional under M1 given the
  individuals already placed (a single Gibbs-like sweep), and recording the
  proposal probability. Dropping IBD indicators for individuals that become
  age 2 is deterministic.
- for each individual sampled in another population C and migrant from A∪B:
  source A or B, drawn from its conditional.

The reverse of a merge is a split that happens to propose the merged state's
preimage, so the Hastings ratio is the product of the sequential proposal
probabilities, divided by the IBD-proposal probability of the merge.

**When A and B are distinct, the move is rejected quickly.** When they are near
identical, the A/B labels are nearly unidentified under M1, and the sequential
conditional proposal is what keeps acceptance high.

## 5. Estimating the Bayes factor

- **Model indicator.** BF_01 = posterior odds of M0 / prior odds.
- **Prior odds tuning.** When the evidence is decisive, the chain rarely visits
  the disfavoured model. Adapt the prior odds during burn-in so both models are
  visited (a pseudo-prior weight), then fix them and correct the posterior
  odds. If one model is never visited, report a bound, as the new BF_01 for m
  does.
- **Rao-Blackwellized estimator.** Average the acceptance probabilities of the
  proposed moves in both directions, in place of counting visits (Bartolucci,
  Scaccia & Mira 2006, for reversible jump). This is more efficient and reuses
  quantities already computed.
- **Partition prior (general form).** Uniform over partitions, or Ewens/CRP
  with a concentration parameter. Results will depend on it, as on any model
  prior, and should be reported with it.

## 6. Caveats

- **Composite likelihood.** Treating linked SNPs as independent overstates the
  information, so the evidence against M0 will be overstated. Thin to one SNP
  per RAD tag (as for the mangrove data), and check calibration by simulation.
- **Exact identity versus practical pooling.** With enough loci, any real pair
  of populations with nonzero differentiation will reject M0. That is the
  correct answer to "are they one deme?". Whether the split model still gives
  stable BA3 estimates for such a pair is a separate question, to be checked by
  simulation (section 7), rather than built into the test.
- **The island model is a model.** Real structure (isolation by distance,
  hierarchical clusters) is represented only through which samples share a
  deme. The test is between island models, not against every structured
  alternative.
- **Interaction with migration.** Under M0, migrants between A and B are
  residents; under M1 near M0, they are unidentifiable. The m prior enters the
  Bayes factor through the Polya-urn terms, as it must, and should be reported.

## 7. Validation

1. **Exact evidence on tiny data.** `exact_posterior.py` and
   `exact_posterior_hier.py` already enumerate the full posterior. log p(G) =
   offset + log Z (check that every prior normalizing constant is included).
   The exact BF for merging A and B is the evidence of the same data with A and
   B relabelled as one population, divided by the evidence of the original.
   Extend `check_exact.py` to compare it with the sampler's estimate.
2. **Python prototype** of the merge/split move in `tests/proto_hier.py`,
   validated against step 1 before any C++.
3. **Simulations** with `sim/simulate_ba3.py` (congen2026 repo):
   - A, B truly one population: simulate one deme and split its sample at
     random into two labels;
   - A, B differentiated at Fst 0.002, 0.005, 0.01, 0.02;
   - A, B with migration to and from a third population;
   - with and without RAD-like linkage (several SNPs per tag), for
     calibration.
   Report Pr(M0), and the m and ancestry estimates under each model.

## 8. Phases

The goal is the full model of section 2: given K population samples, infer how
many demes there are and which samples come from which deme. It is built up
in three steps that share one merge/split move; each step adds only bookkeeping
and a prior.

**Step 1: one fixed pair, K versus K-1 demes.**
Two models only: every sample its own deme (current BA3), or one named pair
(A, B) sharing a deme. The chain merges or splits only that pair.

- 1a. Exact evidence in `exact_posterior.py` and `exact_posterior_hier.py`
  (log p(G) = offset + log Z, with all prior constants), and the exact merge
  Bayes factor: evidence of the data with A and B relabelled as one sample,
  divided by the evidence of the original (section 7.1).
- 1b. Python prototype of the move in `tests/proto_hier.py`: ancestry mapping
  and IBD proposals (section 4), theta_AB by quadrature or proposal, m redrawn
  from its full conditional after every move (partially collapsed Gibbs),
  prior odds adapted during burn-in, and the Rao-Blackwellized Bayes factor
  (section 5). Validate against 1a.
- 1c. C++ with `-c` (and `-H`), e.g. `--pool-test A,B`. Output: Pr(K-1
  demes), BF_01, and estimates conditional on each model. Validate with
  `check_exact.py`.
- 1d. Simulations (section 7.3): a truly oversplit pair, slightly
  differentiated pairs, migration to and from a third deme, thinned and
  linked SNPs.

**Step 2: any one pair, K versus K-1 demes.**
Any pair may merge: K(K-1)/2 merged models plus the split model. Adds the
choice of pair in the move and a prior over pairs (uniform over pairs, with
total prior mass on K-1 demes as one setting). Output: Pr(K-1 demes) and the
posterior probability of each pair being merged. Validate against exact
enumeration over all pairs on tiny data.

**Step 3: full partition model.**
Repeated merges and splits over all partitions of the K samples: merge any two
demes; split a deme into two groups of its samples (proposal over the
bipartitions of its samples). Adds the deme bookkeeping (count tables, IBD
totals, migration rows and theta per deme, sample-to-deme map) and a prior
over partitions (uniform, Ewens/CRP, or uniform on the number of demes then on
partitions). Output: posterior of the number of demes, pairwise co-membership
probabilities, and m, theta and F conditional on the most probable partition
or model-averaged (migration between samples in one deme is undefined, not
zero). Validate against exact enumeration of all partitions of 3-4 samples,
and by simulation.

`--sdpool` (`-P`) and the stepping-stone pooling test on `collapse-m` would be
superseded by this approach for pooling questions. Stepping-stone remains
available for other model comparisons.

## 9. References

- Wright S (1931) Evolution in Mendelian populations. Genetics 16:97-159.
- Balding DJ, Nichols RA (1995) A method for quantifying differentiation
  between populations at multi-allelic loci. Genetica 96:3-12.
- van Dyk DA, Park T (2008) Partially collapsed Gibbs samplers. JASA
  103:790-796.
- Bartolucci F, Scaccia L, Mira A (2006) Efficient Bayes factor estimation from
  the reversible jump output. Biometrika 93:41-52.
- Bryant D, Bouckaert R, Felsenstein J, Rosenberg NA, RoyChoudhury A (2012)
  Inferring species trees directly from biallelic genetic markers. MBE
  29:1917-1932.
- Nicholson G, Smith AV, Jonsson F, et al. (2002) Assessing population
  differentiation and isolation from single-nucleotide polymorphism data.
  JRSS B 64:695-715.
- Pickrell JK, Pritchard JK (2012) Inference of population splits and mixtures
  from genome-wide allele frequency data. PLoS Genet 8:e1002967.
