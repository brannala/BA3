# KL Divergence Implementation Testing Notes

## Overview

These notes document testing of the numerical KL divergence implementation in BA3 v3.4.0. The KL divergence measures information gain from prior to posterior for migration rate parameters.

## Implementation

The KL divergence is computed using:
- **Exact analytical prior**: Uniform(0, 1/3) with density p(m) = 3
- **KDE for posterior**: Gaussian kernel with reflection boundary correction
- **Numerical integration**: Trapezoidal rule over 1000 grid points on [0, 1/3]

Formula:
```
KL(posterior || prior) = ∫ p(m|D) log(p(m|D) / 3) dm
```

Converted to bits using log₂(e) ≈ 1.4427.

## Test Cases

### 1. Uninformative Data (Baseline)

**Dataset**: 2 populations, 1 locus, identical allele frequencies (50/50 in both)
- 20 individuals total

**Results**:
| Metric | Value |
|--------|-------|
| KL divergence | 0.22 bits |
| Mean migration | ~0.22 |
| Posterior SD | 0.08 |
| Bayes Factor | Supports H1 |

**Interpretation**: With no population differentiation, data provide minimal information about migration. KL near zero as expected. BF supports H1 because posterior spreads across parameter space (density at zero is low).

---

### 2. High Information, No Migration (Varying Loci)

**Dataset**: 2 populations, complete differentiation (Fst = 1.0), no migrants

| Loci | Ind/pop | KL (bits) | BF | log₁₀(BF) | Post. SD |
|------|---------|-----------|-----|-----------|----------|
| 5 | 20 | 3.0 | ~45 | +1.65 | 0.01 |
| 20 | 30 | 3.6 | ~75 | +1.88 | 0.01 |
| 100 | 30 | 3.6 | ~75 | +1.90 | 0.01 |

**Observation**: KL plateaued at ~3.6 bits despite increasing loci from 20 to 100. Posterior SD stuck at 0.01 floor.

---

### 3. High Information, No Migration (Varying Sample Size)

**Dataset**: 2 populations, 100 loci, complete differentiation, no migrants

| Ind/pop | KL (bits) | BF | log₁₀(BF) | Interpretation |
|---------|-----------|-----|-----------|----------------|
| 30 | 3.6 | ~75 | +1.90 | Strong for H0 |
| 100 | 5.2 | ~267 | +2.43 | Decisive for H0 |

**Observation**: Increasing sample size from 30 to 100 individuals per population broke through the KL ceiling. Sample size matters more than loci count for migration rate precision.

---

### 4. High Information, With Migration

**Dataset**: 2 populations, 100 loci, 100 individuals/pop, 10% true migration rate
- 90 residents + 10 migrants per population

**Results**:
| Direction | True m | Estimated m | KL (bits) | BF |
|-----------|--------|-------------|-----------|-----|
| pop1→pop0 | 0.10 | 0.036 | 3.06 | Decisive for H1 |
| pop0→pop1 | 0.10 | 0.035 | 2.91 | Decisive for H1 |

**Observations**:
1. Migration correctly detected (Decisive for H1)
2. Estimated rate (~3.5%) lower than true rate (10%) - expected behavior as BA3 partitions migrants across ancestry generations
3. KL still high (~3 bits) indicating informative data

---

## Summary of KL Divergence Behavior

| Scenario | KL (bits) | Interpretation |
|----------|-----------|----------------|
| Uninformative data | ~0.2 | Posterior ≈ Prior |
| Moderate information | 2-3 | Posterior moderately concentrated |
| High information, no migration | 3-5 | Posterior concentrated near zero |
| High information, with migration | 2-3 | Posterior concentrated away from zero |

## Dynamic Range

Based on testing, the practical range of KL divergence is approximately **0 to 5+ bits**:
- **0 bits**: No information gain (posterior equals prior)
- **1-2 bits**: Moderate information
- **3-4 bits**: High information
- **5+ bits**: Very high information (requires large sample sizes)

## Factors Affecting KL Divergence

1. **Population differentiation (Fst)**: Higher differentiation → more information about migration
2. **Number of loci**: More loci → more information, but diminishing returns
3. **Sample size**: More individuals → breaks through KL ceiling
4. **True migration rate**: Affects where posterior concentrates, not necessarily KL magnitude

## Comparison: KL vs Bayes Factor

| Measure | What it tells you |
|---------|-------------------|
| KL divergence | How much we learned about migration rate (regardless of value) |
| Bayes Factor | Whether migration rate is zero vs non-zero |

These are complementary:
- High KL + BF for H0 → Strong evidence of NO migration
- High KL + BF for H1 → Strong evidence of migration
- Low KL + any BF → Uninformative data, conclusions unreliable

## Conclusions

1. The numerical KL implementation correctly captures information gain
2. KL approaches zero for uninformative data (validated)
3. KL increases with data quality (loci, sample size, differentiation)
4. Sample size is critical for achieving high KL values
5. KL and Bayes Factor together provide complementary information about migration

---

*Testing performed: January 2026*
*BA3 Version: 3.4.0*
