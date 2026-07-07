# BA3 sex-biased dispersal — simulation & recovery testing

Ground-truth tooling for the sex-biased dispersal model
(`doc/sex_biased_dispersal_model.md`). Pure Python 3 stdlib, no dependencies.

## Workflow

```
# 1. Simulate data. --phi sets both sex ratios; use --phi-move / --phi-breed to
#    make movement (age-1) and gene flow (age-2) differ, e.g. males move but
#    females' gene flow dominates:
python3 sim/simulate_sexbias.py --phi-move 0.3 --phi-breed 0.7 --m 0.12 \
    --nauto 200 --nx 250 --nind 250 --fst 0.15 --seed 5 --out sim_run

# 2. Run BA3 (VCF + metadata; sex-biased model turns on when sex is present)
./BA3 -V sim_run.vcf -M sim_run_meta.txt -o BA3out.txt

# 3. Score BA3's estimates (phiMove, phiBreed, rho) against the ground truth
python3 sim/score_recovery.py --truth sim_run_truth.json --ba3out BA3out.txt
```

`score_recovery.py` exits 0 (PASS) if every compared quantity is within `--tol`
standard deviations of its target, else 1 — so it can gate an automated test.

## Files produced by the simulator

| file | contents |
|---|---|
| `<out>.vcf` | `chr1` autosomal + `chrX` X-linked biallelic SNPs; males haploid on `chrX` |
| `<out>_meta.txt` | `sampleID <tab> pop <tab> sex` (sex = `F`/`M`) |
| `<out>_truth.json` | true `phi_move`/`phi_breed`, migration rate, allele freqs, per-individual ancestry |

## Output contract (parsed by `score_recovery.py`)

BA3 emits a sex-biased dispersal block whose parameter lines begin with the
parameter name and carry a `NUMBER(SD)` token, e.g.

```
   rho      residents (non-dispersers)     = 0.4889(0.0239)
   phiMove  age-1 migrants (movement)      = 0.2248(0.0525)   P(>rho) = 0.000
   phiBreed age-2 migrant parents (gene flow) = 0.6620(0.0576)   P(>rho) = 0.996
```

The scorer compares `phiMove`/`phiBreed` to the simulated truth and `rho` to the
realized resident ratio, and exits 0 (PASS) if each is within `--tol` SDs.

## Notes

- The simulator uses `migrantAge in {0,1,2}` exactly as BA3 does: non-migrant,
  first-generation migrant, second-generation (F1 hybrid) migrant. No
  third-generation state, so no recombinant X mosaic (see the design doc).
- **Movement vs gene flow:** `phiMove` (age-1 migrants) is who *moved*;
  `phiBreed` (age-2 migrant parents) is who moved *and bred* = effective gene
  flow. The X markers inform `phiBreed` (via age-2 males). For tight recovery,
  raise `--m`, `--nind`, and `--nx` so there are enough migrants of each class.
- v1 is HWE within populations (`F = 0`) and biallelic SNPs only.
