# BA3 Project Context

## What is BA3?

BA3 (BayesAss Edition 3) is a Bayesian population genetics program that estimates recent migration rates between populations using multilocus genotype data. It uses MCMC sampling to infer:

- Migration rate matrices between populations
- First and second generation migrant individuals
- Population allele frequencies
- Inbreeding coefficients

## Quick Reference

### Build
```bash
make clean && make
```

### Run
```bash
./BA3 [options] <input_file>
./BA3 -v -M metadata.txt input.vcf  # For VCF files
```

### Key Files
- `src/main.cpp` - Main implementation
- `include/BA3.h` - Data structures and declarations
- `examples/` - Test datasets

### Dependencies
- GSL (GNU Scientific Library)
- htslib (for VCF support)

## Code Conventions

- C++11 standard
- Global variables for MCMC state
- Functions prefixed by purpose: `propose*`, `update*`, `compute*`, `read*`
- GSL random number generator used throughout

## When Modifying Code

1. Test with example datasets after changes
2. Verify MCMC acceptance rates remain reasonable (20-60%)
3. Check output format compatibility if changing output functions
4. Update version string in main.cpp for releases
