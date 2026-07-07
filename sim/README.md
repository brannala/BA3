# BA3 sex-biased dispersal — simulation & recovery testing

Ground-truth tooling for the sex-biased dispersal model
(`doc/sex_biased_dispersal_model.md`). Pure Python 3 stdlib, no dependencies.

## Workflow

```
# 1. Simulate data with known phi and migration rates
python3 sim/simulate_sexbias.py --phi 0.8 --m1 0.10 --m2 0.05 \
    --nauto 300 --nx 150 --nind 60 --fst 0.1 --seed 42 --out sim_run

# 2. Run BA3 on the simulated VCF + metadata (once X/sex support is built)
./BA3 -v -M sim_run_meta.txt sim_run.vcf

# 3. Score BA3's estimates against the ground truth
python3 sim/score_recovery.py --truth sim_run_truth.json --ba3out BA3out.txt
```

`score_recovery.py` exits 0 (PASS) if every compared quantity is within `--tol`
standard deviations of its target, else 1 — so it can gate an automated test.

## Files produced by the simulator

| file | contents |
|---|---|
| `<out>.vcf` | `chr1` autosomal + `chrX` X-linked biallelic SNPs; males haploid on `chrX` |
| `<out>_meta.txt` | `sampleID <tab> pop <tab> sex` (sex = `F`/`M`) |
| `<out>_truth.json` | true `phi`, migration rates, allele freqs, and per-individual true ancestry |

## Output contract for `phi`

`score_recovery.py` looks for a line in `BA3out.txt` mentioning `phi` with a
`NUMBER(SD)` token, e.g.

```
 Dispersal female fraction phi: 0.7912(0.0345)
```

When the BA3 sex-bias code is added, emit `phi` in that form (any line containing
`phi` and a `mean(sd)` token is accepted).

## Notes

- The simulator uses `migrantAge in {0,1,2}` exactly as BA3 does: non-migrant,
  first-generation migrant, second-generation (F1 hybrid) migrant. No
  third-generation state, so no recombinant X mosaic (see the design doc).
- `phi` is informed by the sex of every first-generation migrant: the age-1
  individuals' own (observed) sex, plus the age-2 individuals' migrant-parent
  sex, with the X markers adding the age-2 **male** signal. For tight `phi`
  recovery, raise `--m1/--m2`, `--nind`, and `--nx` so there are enough migrants
  — a handful of lineages gives a noisy empirical `phi`.
- v1 is HWE within populations (`F = 0`) and biallelic SNPs only.
