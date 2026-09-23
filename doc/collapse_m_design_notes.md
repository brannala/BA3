# Collapsed samplers, model comparison and linkage: design notes

Status as of 2026-09-23. Branch `collapse-m` (on top of the 3.5.0 release).

## 1. Summary

- **Released (3.5.0, `master`)**: the collapsed sampler `-c` (allele frequencies
  integrated out), and fixes to two long-standing errors in the default sampler
  (missing-genotype imputation; allele-frequency move for loci with 3+ alleles).
  Validated against exact posteriors.
- **Experimental (`collapse-m`, not released)**: `-C` also integrates out the
  migration rates, with Gibbs ancestry updates and Rao-Blackwellized output;
  `-S` stepping-stone marginal likelihood; `-P` Savage-Dickey pooling test. All
  validated against exact values on small data.
- **Main open problem**: on large SNP data (e.g. `litt`, 12,718 SNPs) the
  likelihood assumes linkage equilibrium among all SNPs. It is a composite
  likelihood that greatly overstates the information in the data, which makes
  posteriors overconfident, creates local modes, and makes tempering and
  stepping-stone very expensive. The planned fix is a **haplotype-block model**
  (section 5).

## 2. Released in 3.5.0

| Item | Notes |
|---|---|
| `-c/--collapse` | Dirichlet-multinomial marginal over allele frequencies; count tables, O(1) ancestry-category selection, incremental migrant counts, log-integer table, memoized homozygote IBD log. F by a Gibbs step on IBD indicators: one individual's indicators refreshed per iteration, F drawn from its Beta full-conditional every 10 iterations, start from a moment estimate. Missing genotypes marginalized exactly. |
| Fix: missing-genotype MH | Ordered-pair proposal was scored with the unordered heterozygote probability (factor 2), double-weighting imputed heterozygotes. Present in every version in the repository. |
| Fix: allele-frequency move | Missing Jacobian (K-2) log[(1-p')/(1-p)] for K >= 3 alleles; proposals above 1 wrapped instead of reflected; possible `nan` frequencies. Present in every version in the repository (boundary code commented 2011). The analogous F-move wrap bug was fixed in 3.0.4 (2018). |
| Tests | `tests/regress_default.sh` (byte identity of the default path), `check_exact.py` / `exact_posterior.py` (exact posterior by enumeration), `compare_collapse.py`, `bench.py`. |

Distribution: GitHub release with binaries; Homebrew tap updated; bioconda PR
bioconda/bioconda-recipes#69544 passed CI and awaits a reviewer.

## 3. Experimental work on `collapse-m`

| Option | Design | Validation |
|---|---|---|
| `-C/--collapse-m` | With u = 3m the uniform prior on each migration row is a flat Dirichlet, so m integrates out given ancestry counts. Each iteration draws one individual's ancestry exactly from its 2P-1 states: Polya-urn prior (n0+1, (c_j+1)/3, 2(c_j+1)/3) times the collapsed genotype predictive. m mean/SD Rao-Blackwellized; ancestry probabilities averaged from stored full-conditionals. Trace and kernel Savage-Dickey use draws of m. | `check_exact.py` passes (m, F, frequencies, ancestry, means and SDs). |
| `-S/--ss` | Stepping-stone, 16 rungs at (k/16)^(1/0.3), 40% burn-in per rung. Parameters: ancestry, IBD indicators of homozygotes, age-2 phase; priors at power 1. Powered likelihood = Dirichlet-multinomial of gene counts x (1-F) per age-0/1 heterozygote (the heterozygote factor must be in the likelihood, otherwise the beta = 0 rung does not sample the prior: this caused a 1.6-2.5 nat bias in the first version). Reports per-rung mean and SD of the log-likelihood. | `check_evidence.py`: log p(G) within 0.005 nats of exact on six cases; pooled-vs-full log BF +0.431 +- 0.004 (exact +0.433). Small systematic bias ~ +0.003 nats. |
| `-P/--sdpool` | Savage-Dickey BF01 for H0: equal allele frequencies of a population pair (nested; not the pooled model), Rao-Blackwellized via the closed-form overlap of Dirichlet full-conditionals. Reports SD of the log term as a reliability diagnostic. | log10 BF01 within 0.002 of exact on three pairs. Estimator is a posterior mean of a product over loci: unstable when SD(ln ratio) is large, i.e. with many loci. |

Performance: on `3pop`, `-C` costs about the same per iteration as `-c` and
gives 3.5x more effective samples per second for the slowest migration rate.
On `litt` (P = 14), `-C` is 2.7x slower per iteration, but for 175 of 182
migration rates its Rao-Blackwellized estimates are identical to 4 decimals
across seeds after 200k iterations, a lower bound of ~700x more effective
samples per second than `-c`.

## 4. Findings on `litt` (106 individuals, 14 populations, 12,718 SNPs)

- **Poorly mixing cluster**: seven rates among pop10, pop11, pop14 (BA3 indices
  9, 10, 13) and pop3 differ between seeds in both `-c` and `-C`; seeds settle in
  different configurations (e.g. mostly pop14<-pop10 vs pop10<-pop14).
- **Not near-identical populations**: Hudson F_ST pop10-pop11 0.22,
  pop11-pop14 0.24, pop10-pop14 0.59 (median over all pairs 0.50). pop11 is
  intermediate between pop10 and pop14. Sample sizes are 8 each (range 5-9).
- **Interpretation**: with ~8 individuals per population and m ~ 0.1, a mode
  corresponds to 2-3 individuals per population assigned as migrants or
  residents; because frequencies are defined by who is assigned, and 12,718
  SNPs make each individual's assignment near-deterministic given the others,
  single-individual moves cannot cross between configurations.
- **Flip diagnostic**: 6 seeds of `-C -g` were run on loki
  (`~/ba3_fixcheck/flip/`); analysis pending.
- **Stepping-stone pooling test (pooled pairs/triple vs full)**: inconclusive.
  Seed-to-seed SD of log p(G) up to ~1,200 nats, larger than the differences.
- **Tempering ladder** (from per-rung statistics, 2 runs): SD(logL) ~ 80 for
  beta in [0.06, 1] (well equilibrated), but 18,000-28,000 for beta < 0.001;
  mean logL rises from ~ -830,000 (prior) to ~ -300,000, almost all below
  beta = 0.004. Implied MC^3 ladder: ~51 chains to beta_min = 0.1 and ~59 to
  0.01 at 30% swap acceptance; communication barrier Lambda ~ 66, i.e. ~130
  chains for non-reversible PT from the prior. This also explains the
  stepping-stone failure: the fixed schedule has only 3 rungs below 0.004,
  with dbeta * SD ~ 50 on one step.

## 5. Diagnosis: linkage and the composite likelihood

- The 12,718 SNPs lie on 5,201 RAD tags (mean 2.45 SNPs/tag, up to 16; 84% of
  SNPs share a tag). All loci have CHROM = `chr1`, so there are no positions.
- For **unphased** data, BA3's migrant categories are exact under linkage
  equilibrium within populations (a migrant carries two source haplotypes; an
  age-2 individual one source and one native haplotype, and its per-locus phase
  sum is then the exact chromosome-level likelihood). The error is the
  assumption of **linkage equilibrium within populations**: linked SNPs carry
  largely redundant information and the product over loci counts it
  repeatedly. Thinning would reduce the peakedness but is not the solution.
- Consequences: overconfident ancestry and migration posteriors, deep local
  modes, large SD(logL) and Lambda, unstable stepping-stone estimates, and
  inflated Bayes factors for population partitions.

## 6. Plans

### 6.1 Haplotype-block model (priority)

- Loci become **blocks** with population-specific **haplotype frequencies**
  (arbitrary LD within a block, blocks independent). Resident: two native
  haplotypes; first-generation migrant: two source haplotypes; age 2: one of
  each (origin latent per block; shared along a chromosome if positions are
  available). IBD/F at the block level.
- **Collapsed**: Dirichlet (or Dirichlet-process) prior on haplotype
  frequencies, integrated out; a block is a multi-allelic locus whose alleles
  are haplotypes. Sparse counts (hash of observed haplotypes per population and
  block). Existing machinery (urn ancestry Gibbs, m integration, stepping-stone,
  exact tests) generalizes.
- **Phase latent**: each individual's pair of haplotypes per block (2^(h-1)
  pairs for h heterozygous sites) Gibbs-sampled given the others' counts and
  summed inside the ancestry update; missing SNPs within a block are latent
  alleles.
- **Open decisions**:
  1. Blocks: RAD tags from the ID prefix (5,201 blocks for `litt`), or windows
     along chromosomes when positions exist (long-range LD would then need a
     Li-Stephens / fastPHASE-type mosaic model: a larger project).
  2. Prior: symmetric Dirichlet over 2^k haplotypes with concentration theta, or
     a Dirichlet process with a linkage-equilibrium base measure (nests the
     current model as theta -> infinity; base allele frequencies fixed at pooled
     estimates or given their own prior). theta fixed with sensitivity analysis,
     or estimated.
- **Validation**: extend `exact_posterior.py` to small blocks (k <= 3, summing
  over phase); `check_exact.py` / `check_evidence.py`; then rerun the `litt`
  ladder estimate as a benchmark (SD(logL) and Lambda should fall).

### 6.2 Mixing across configurations

- **Adaptive non-reversible parallel tempering** (Syed et al. 2022): beta = 0
  reference chain sampled exactly from the prior (urn ancestry, F, IBD); adapt
  the schedule during burn-in to equalize swap rejection and set N ~ 2 Lambda
  (user cap); chains in threads, per-chain RNG, swaps of beta labels at
  barriers; output from the beta = 1 chain. Requires moving chain state into a
  struct (count tables, IBD totals, ancestry, F, RNG, ssBeta).
- **Group reassignment move**: e.g. switch all migrants from source j into s to
  source k, or move a whole (s, j, age) category; exact MH with the urn prior
  and collapsed likelihood. Cheaper than tempering for the specific modes seen;
  combines with it.
- Decide after the haplotype model, since Lambda under the current likelihood
  (~66) is probably inflated by the linkage problem.

### 6.3 Model comparison

- **Adaptive stepping-stone schedule**: place rungs so dbeta * SD(logL) is
  about constant (pilot run as in section 4); shares code with NRPT.
- **Rao-Blackwellized Savage-Dickey for m = 0** under `-C`: exact density at 0
  is 1[c_j = 0] * 3(N_s + P - 1) against prior 3(P - 1); replaces the kernel
  estimate but changes reported BF/KL values (decision pending).
- Savage-Dickey pooling (`-P`) remains a screen; stepping-stone is the formal
  test of pooled vs separate populations.

### 6.4 Housekeeping

- Merge `collapse-m` into `master` once the haplotype model and mixing decisions
  settle; document `-C`, `-S`, `-P` in the manual and wiki.
- `compare_collapse.py`: ESS-based standard errors (between-seed SEs are
  optimistic for slowly mixing quantities).
- Benchmarks still to run: N-scaling of `-c`/`-C`; `vcf_litt` against a long
  standard-sampler reference.
- `dev` still carries the sex-biased dispersal / gamma features (not in
  `master`); decide how they relate to the collapsed samplers.

## 7. Pointers

- Code: `src/main.cpp` (`computeAddLogProb[Beta]`, `addIndividual`,
  `resampleIBD`, `drawF`, `migRowMoments`, `ssInit`/`ssAccum`,
  `sdPoolLogRatio`); tests in `tests/` (see `tests/README.md`).
- loki: `~/ba3_fixcheck/` (builds and runs; `littss/` stepping-stone model
  comparison, `ladder/` per-rung runs, `flip/` flip diagnostic).
