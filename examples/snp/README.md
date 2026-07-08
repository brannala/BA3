# Large empirical SNP example datasets

Published, openly-licensed (CC0) multi-population SNP datasets for testing BA3 —
especially the collapsed (`--collapse`) sampler at realistic locus counts. BA3
reads the VCFs natively (`-V vcf -M meta`); no format conversion is needed.

The VCFs and their `INDIV POPLN` metadata are included here so the examples run
out of the box. `./fetch_snp_examples.sh` re-downloads them from Dryad and
regenerates the metadata (useful to refresh or to see exactly how they were
prepared).

## Datasets

### `dace/` — Speckled dace, Greater Death Valley (*Rhinichthys osculus*)
- **145 individuals, 13 populations, 14,355 SNPs** (`spd_dv.filt2.vcf`), ~28% missing.
- ddRAD; population = the letter code in each sample name (e.g. `1AMA001` → `AMA`,
  `1damr001` → `DAMR`).
- Dryad: <https://doi.org/10.5061/dryad.51c59zw62> (file `spd_dv_input_files.tar.gz`).
- Chafin, T. K., Douglas, M. R., Martin, B. T., & Douglas, M. E. (2019). Hybridization
  drives genetic erosion in sympatric desert fishes of western North America — see the
  Death Valley speckled dace study, *Heredity*/*Mol. Ecol.* series (Dryad above).

### `palm/` — Malagasy palms (4 species)
Per-species VCFs with missing data retained (`*_miss50.vcf`, up to 50% missing).
Population = the sample-name prefix before the first `_` (e.g. `HC1_LM007` → `HC1`).

| Species | File | Indiv | Pops | SNPs |
|---|---|---|---|---|
| *Bm* | `Bm_miss50.vcf` | 27 | 4 | 1,935 |
| *BN* | `BN_miss50.vcf` | 52 | 8 | 5,524 |
| *Cm* | `Cm_miss50.vcf` | 45 | 7 | 11,602 |
| *HC* | `HC_miss50.vcf` | 43 | 7 | 12,913 |

- Dryad: <https://doi.org/10.5061/dryad.pc866t1z6> (megafrugivore-dispersal palms).
  The archive also ships the authors' own BayesAss `.immanc` inputs under
  `input/BayesAss/`.

## Running

Base migration model (these have no sex/X information):

```bash
# standard sampler
./BA3 -V examples/snp/dace/spd_dv.filt2.vcf -M examples/snp/dace/spd_dv_meta.txt \
      -i 5000000 -b 1000000 -n 2000 -o dace_std.txt

# collapsed sampler (integrates out allele frequencies; faster + better mixing)
./BA3 -c -V examples/snp/dace/spd_dv.filt2.vcf -M examples/snp/dace/spd_dv_meta.txt \
      -i 5000000 -b 1000000 -n 2000 -o dace_col.txt
```

Notes
- Missing genotypes: the standard sampler imputes them; `--collapse` marginalizes
  them out (drops the missing gene copies). The palm `miss50` sets (and the ~28%-
  missing dace) exercise that path — expect the two samplers to differ somewhat on
  heavily-missing data, with collapse the less biased (see
  `doc/collapsed_allele_freqs_plan.md`).
- These are ddRAD SNPs on RAD-tag contigs (no `chrX`), so the sex-biased model is
  not engaged; they test the base migration model at scale.
- For a mixing/effective-sample-size comparison, see `sim/run_mixing.sh`.
