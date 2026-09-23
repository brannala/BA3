# BA3 tests

Validation and benchmark scripts for the collapsed sampler (`-c/--collapse`) and
for the standard sampler. The Python scripts need numpy (`exact_posterior.py`,
`check_exact.py`, `bench.py`).

```bash
# a reference build, e.g. the last release
git worktree add /tmp/ba3-ref v3.5.0 && make -C /tmp/ba3-ref

tests/regress_default.sh /tmp/ba3-ref/BA3 ./BA3              # ~4 min
tests/check_exact.py ref=/tmp/ba3-ref/BA3 collapse=./BA3:-c  # exact posterior; ~5 min on 8 cores
tests/compare_collapse.py /tmp/ba3-ref/BA3 ./BA3              # ~1 h on 8 cores (--quick: ~5 min)
tests/bench.py /tmp/ba3-ref/BA3 ./BA3                         # run on an idle machine
```

| Script | Checks |
|---|---|
| `regress_default.sh` | Without `-c`, the new binary is **byte-identical** to the reference for the same seed: all output files (main output, trace, `BA3indiv.txt`, `-F` frequencies) and stdout, on native and VCF inputs, with and without missing data. |
| `check_exact.py` | Every build given is run with R seeds on the tiny data sets `data/exact_*.txt`, and each posterior mean and SD (migration rates, F, allele frequencies, individual ancestry probabilities) is compared with the **exact posterior** from `exact_posterior.py`. A group fails if any \|z\| > 4 or more than 2% exceed 3. |
| `exact_posterior.py` | Exact BA3 posterior by enumeration of all (2P-1)^N ancestry assignments, with migration rates (Dirichlet), allele frequencies (Dirichlet-multinomial, IBD indicators and age-2 phase summed) and F (Gauss-Legendre quadrature, exact for these polynomial integrands) integrated out. Keep N small (8 individuals, 2 populations: 6561 assignments, ~4 s). |
| `check_evidence.py` | Model-comparison estimators against exact values: the stepping-stone log marginal likelihood (`-S`) on every exact data set and on `exact_3pop` with two populations pooled, the pooled-vs-full log Bayes factor, and the Savage-Dickey pooling test (`-P`, H0: equal allele frequencies of a population pair). Exact values come from `exact_posterior.py --pool` / `--tie`. |
| `compare_collapse.py` | With `-c`, posterior means agree with a standard-sampler build on the example data sets (R seeds per mode, z-scores from between-seed standard errors). Between-seed SEs are optimistic for slowly mixing quantities, so use long chains. |
| `bench.py` | Seconds per 100k iterations, speedup, effective sample size (ESS) of the migration-rate traces, and ESS per second for both samplers. `--snp DIR` adds the empirical SNP sets. |

The exact data sets each target one feature (`exact_3pop`, with three populations, is for the pooling tests): `exact_k2` (biallelic loci),
`exact_k234` (loci with 2, 3 and 4 alleles), `exact_missing` (`exact_k234` with
three missing genotypes) and `exact_highF` (all homozygous).

**Known reference problems.** The standard sampler in BA3 up to 3.4.4 had two
bugs, fixed in 3.5.0, that `check_exact.py` detects: the missing-genotype update over-weighted
imputed heterozygotes by a factor of 2, and the allele-frequency move lacked the
Jacobian for loci with 3 or more alleles (and wrapped rather than reflected at 1,
which could also produce `nan` frequencies). A 3.4.x build therefore
fails `check_exact.py` on `exact_k234`, `exact_missing` and `exact_highF`, and
as the reference for `compare_collapse.py` it produces disagreements on data with
missing genotypes or multi-allele loci that come from the reference, not from
`-c`. `regress_default.sh` is for comparing against a build with the same
standard sampler (3.5.0 or later): the reflection fix changes the chain, though
not the posterior, even for biallelic loci.

`data/2pop_missing.txt` and `data/100loci_missing.txt` are missing-data variants
of example files, made with `make_missing.py` (10% of genotypes fully missing,
2% half-missing, seed 1).

Running `BA3 -c -d` also prints two self-checks: removing and re-adding each
individual's gene copies incrementally must reproduce the full recomputed
Dirichlet-multinomial log-likelihood, and the running IBD/genotype totals used
for the inbreeding update must match a recount at the end of the run.
