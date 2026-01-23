# BA3 (BayesAss Edition 3) Development Skill

## Project Overview

BA3 is a Bayesian population genetics tool for estimating recent migration rates between populations using multilocus genotype data. Developed by Bruce Rannala at UC Davis, it implements Markov Chain Monte Carlo (MCMC) inference to estimate:

- Per-generation migration rates between populations
- Individual migrant ancestry (first and second generation migrants)
- Population allele frequencies
- Inbreeding coefficients (F-statistics)

**Version**: 3.4.0 (December 2025)
**License**: GNU Affero General Public License v3
**Language**: C++11

## Codebase Structure

```
BA3/
├── src/main.cpp          # Core implementation (~3,000 lines)
├── include/BA3.h         # Header with data structures and declarations
├── examples/             # Example input files
│   ├── 2pop.txt          # 2-population example
│   ├── 3pop.txt          # 3-population example
│   └── allpopslitt_meta.txt  # VCF metadata example
├── Makefile              # Build configuration
└── README.md             # User documentation
```

## Key Data Structures (include/BA3.h)

- `struct indiv`: Individual genotype data, migrant ancestry, likelihood
- `struct ancestryProbs`: Posterior probabilities for ancestry categories
- `struct SavageDickeyStats`: Statistics for Bayesian hypothesis testing

## Dependencies

| Library | Purpose | Installation |
|---------|---------|--------------|
| GSL | Random numbers, gamma functions | `apt install libgsl-dev` or `brew install gsl` |
| htslib | VCF file support | `apt install libhts-dev` or `brew install htslib` |

## Building

```bash
# Standard build
make clean && make

# Debug build
make clean && make DEBUG=1

# Static build (Linux)
make clean && make STATIC=1
```

## Theoretical Foundation

### Bayesian Model

BA3 estimates the posterior distribution:
```
P(m, p, F | G) ∝ P(G | m, p, F) × P(m) × P(p) × P(F)
```

Where:
- `m[i][j]` = migration rate from population j to i
- `p` = allele frequencies per population/locus
- `F` = inbreeding coefficients per population
- `G` = observed genotype data

### MCMC Proposals

1. **Migration rates**: Random walk with reflection at boundaries
2. **Individual ancestry**: Add/drop migrants with Metropolis-Hastings
3. **Allele frequencies**: Log-space random walk
4. **Inbreeding coefficients**: Random walk on [0,1]
5. **Missing genotypes**: Sample from conditional distribution

### Constraints

- Diagonal migration rates: `m[i][i] >= 2/3` (most individuals are non-migrants)
- Off-diagonal sum: `Σ m[i][j] <= 1/3` for j ≠ i

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `-i N` | 500,000 | Total MCMC iterations |
| `-b N` | 250,000 | Burn-in iterations |
| `-n N` | 100 | Sampling interval |
| `-m D` | 0.10 | Migration rate proposal step |
| `-a D` | 0.10 | Allele frequency proposal step |
| `-f D` | 0.10 | F-statistic proposal step |
| `-s N` | 10 | Random seed |
| `-t` | off | Output trace file |
| `-g` | off | Output individual ancestry |
| `-v` | off | Read VCF format |
| `-o FILE` | BA3out.txt | Output file name |

## Input Formats

### BA3 Native Format
```
<n_loci> <n_alleles_1> <n_alleles_2> ...
<ind_id> <pop> <locus1_a1> <locus1_a2> <locus2_a1> <locus2_a2> ...
```

### VCF Format (with -v flag)
Requires metadata file (-M) mapping samples to populations:
```
sample_id<TAB>population_name
```

## Output Files

1. **BA3out.txt**: Migration matrices, F-statistics, allele frequencies
2. **BA3indiv.txt** (-g): Individual ancestry posterior probabilities
3. **BA3trace.txt** (-t): MCMC trace for convergence diagnostics

## Common Development Tasks

### Adding a new command-line option
1. Add variable declaration in main.cpp globals
2. Add to getopt_long options array
3. Add case in switch statement
4. Update help text in usage()

### Modifying MCMC proposals
- Migration: Look for `proposeM()` function
- Ancestry: Look for `proposeAncestry()` function
- Allele frequencies: Look for `proposeA()` function

### Adding new output
1. Modify output functions at end of main()
2. Update header file if new data structures needed

## Convergence Diagnostics

Good MCMC runs show:
- Acceptance rates: 20-60% for each proposal type
- Trace plots: Stable plateau after burn-in
- Multiple runs with different seeds: Consistent posteriors

Auto-tuning during burn-in adjusts delta values to achieve target acceptance rates.

## Savage-Dickey Tests (v3.4.0)

Tests H₀: m[i][j] = 0 (no migration) vs H₁: m[i][j] > 0

Output includes:
- Bayes factors with Jeffreys scale interpretation
- Kullback-Leibler divergence
- Prior and posterior density estimates at zero

## Testing

```bash
# Run with example data
./BA3 -i 10000 -b 5000 -n 10 examples/3pop.txt

# Check convergence with trace output
./BA3 -t -i 100000 -b 50000 examples/3pop.txt
```

## Debugging Tips

1. Use `-s` option to set reproducible seed
2. Enable trace output (`-t`) to monitor MCMC
3. Check acceptance rates in terminal output
4. For VCF issues, verify metadata file format matches sample IDs

## Performance Considerations

- Likelihood calculation is O(individuals × loci)
- Memory scales with: populations × loci × max_alleles
- Large datasets: Consider longer burn-in, larger sampling intervals
