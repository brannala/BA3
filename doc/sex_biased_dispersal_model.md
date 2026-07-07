# Sex-Biased Dispersal Model for BA3

**Status:** design draft
**Scope:** extends BA3's migrant-ancestry model to infer a global dispersal
sex-bias from autosomal + X-linked markers with per-individual sex metadata,
keeping BA3's existing two-generation ancestry depth (`migrantAge` 0, 1, 2).

---

## 1. Goal

BA3 currently estimates a single directional migration rate `M_{s<-b}` (fraction
of population `s` that is migrant-derived from source `b` per generation). This
extension splits each rate by the **sex of the dispersing (first-generation)
migrant** and infers a single global **dispersal sex-bias**, while retaining all
standard BA3 assumptions (HWE within populations, constant allele frequencies,
unlinked loci, inbreeding coefficient `F` per population) and BA3's existing
ancestry depth.

The signal comes from **X-linked inheritance being sex-asymmetric**: a male
transmits his single X only to daughters, and a second-generation individual's X
therefore records whether its migrant parent was the mother or the father. This
lets the data inform the sex composition of migrants.

---

## 2. BA3 ancestry states (unchanged)

Each individual carries a latent `migrantAge` (BA3's existing state), traced back
at most two generations:

| `migrantAge` | Meaning | Autosomal genotype model |
|---|---|---|
| **0** | non-migrant (native resident) | both copies from the native population `s` (with `F_s`) |
| **1** | first-generation migrant (individual migrated) | both copies from source `b` (with `F_b`) |
| **2** | second-generation migrant (a parent migrated; F1 hybrid) | one copy from `b`, one from `s`; phase summed |

There is no third-generation (grandparent-migrant) state. A consequence used
throughout: **no X-chromosome recombinant mosaic ever arises**, because a mosaic
X would require the focal individual's mother to herself be a second-generation
migrant (carrying one source X and one native X) — which would make the focal a
third-generation migrant, outside BA3's depth. At `migrantAge <= 2` every X copy
is a whole chromosome, either all-source or all-native.

---

## 3. Parametrization (female-fraction form)

Split each baseline rate `M_{s<-b}` by the sex of the first-generation migrant
using **one global female fraction** `phi in (0,1)`:

```
M_{s<-b}^female = M_{s<-b} * phi
M_{s<-b}^male   = M_{s<-b} * (1 - phi)
```

Interpretation:

- `phi` = fraction of migrants that are female (global, shared across all pairs
  and directions).
- `phi > 1/2` female-biased dispersal, `phi < 1/2` male-biased, `phi = 1/2` no bias.
- The sex-summed rate is preserved: `M^female + M^male = M_{s<-b}`, so **every
  `M_{s<-b}` keeps its current meaning**, BA3's feasibility constraint
  `sum_{b != s} M_{s<-b} <= 1` is unchanged, and existing priors/output formats
  carry over. (This is why the fraction is preferred over a multiplicative ratio,
  which would inflate the row sums and couple the bias into the constraint.)
- Bias ratio for reporting: `r = phi / (1 - phi)`.

**Parameter count:** for `P` populations, `2 P(P-1)` sex-specific rates collapse
to `P(P-1) + 1` (the `M`'s plus `phi`). For 2 populations: 4 -> 3. Testing
`phi = 1/2` (no sex bias) is a Savage-Dickey density-ratio test, directly
analogous to BA3's existing `M = 0` tests.

---

## 4. Data and assumptions

Per individual `i`:

- `s_i`  : sampled population (observed).
- `g_i`  : sex, `g_i in {F, M}` (observed metadata).
- Genotypes at autosomal loci `L_A`: diploid ordered pair `(u, v)`.
- Genotypes at X-linked loci `L_X`:
  - `g_i = F`: diploid `(u, v)`.
  - `g_i = M`: hemizygous, single allele `u`.

Loci are annotated as autosomal or X-linked; **pseudoautosomal regions are
excluded** from `L_X`. Standard BA3 assumptions retained: within-population HWE
with inbreeding `F_p`, constant allele frequencies `f_{p,l,a}`, loci unlinked,
and migrant ancestry traced back at most two generations.

---

## 5. Parameters

| Symbol | Meaning |
|---|---|
| `f_{p,l,a}` | frequency of allele `a` at locus `l` in population `p` (autosomal and X) |
| `F_p` | inbreeding coefficient of population `p` |
| `M_{s<-b}` | baseline (sex-summed) migration rate into `s` from `b` |
| `phi` | global female fraction of migrants |

---

## 6. Latent state

Each individual carries `(b, migrantAge, sigma)`:

- `b`          : source population of the migrant lineage (`b != s_i`); absent for `migrantAge = 0`.
- `migrantAge` : `0`, `1`, or `2` (Section 2).
- `sigma`      : the **new** label — the sex of the first-generation migrant.
  - `migrantAge = 0`: none.
  - `migrantAge = 1`: `sigma = g_i` (the migrant is the focal individual; its sex
    is observed, so `sigma` is fixed, not free).
  - `migrantAge = 2`: `sigma in {F, M}` = the sex of the migrant parent (latent).

This adds a single binary label to at most the `migrantAge = 2` individuals.

---

## 7. Genotype likelihood building blocks

Let `s = s_i` (native background) and `b` the migrant source.

Within-population diploid genotype probability with inbreeding:

```
H_p(u,v) = F_p * f_{p,u} * 1[u=v] + (1 - F_p) * f_{p,u} * f_{p,v}
```

Single X gene-copy allele probability, for a copy of type `tau` observing allele
`a` at X-locus `l`:

```
psi_l(a | B) = f_{b,l,a}     (all-source copy)
psi_l(a | N) = f_{s,l,a}     (all-native copy)
```

(There is no mosaic type at this ancestry depth; see Section 2.)

### 7.1 Autosomal loci (unchanged from BA3)

Per autosomal locus:

```
migrantAge = 0 : H_s(u,v)
migrantAge = 1 : H_b(u,v)
migrantAge = 2 : (1/2)[ f_{b,u} f_{s,v} + f_{b,v} f_{s,u} ]   (outbred, no F)
```

This is exactly the current `logLik` (`main.cpp` age 0 / 1 / 2 branches).

### 7.2 X-linked loci

X copy-types by state:

| `migrantAge` | focal sex | condition | X copy-type(s) | reads migrant sex? |
|---|---|---|---|---|
| 0 | M | — | `{N}` | — |
| 0 | F | — | `{N,N}` (apply `F_s`) | — |
| 1 | M | — | `{B}` | migrant is focal (observed) |
| 1 | F | — | `{B,B}` (apply `F_b`) | migrant is focal (observed) |
| 2 | M | mother migrated (`sigma = F`) | `{B}` (whole X) | **yes** |
| 2 | M | father migrated (`sigma = M`) | `{N}` (whole X) | **yes** |
| 2 | F | either `sigma` | `{B,N}` | no (identical both ways) |

Rationale for the `migrantAge = 2` rows:

- A male gets his single X from his mother. If the mother is the migrant she is
  entirely source, so his X is all-source `{B}`; if the father is the migrant the
  father contributes only the Y and the X comes from the native mother, so his X
  is all-native `{N}`. Hence a **2nd-gen male's X directly reads the migrant
  parent's sex.**
- A female gets one X from each parent. Whichever parent is the migrant
  contributes an all-source X and the native parent an all-native X, so she is
  `{B,N}` **regardless of which parent migrated** — her X carries no
  sex-of-migrant information (though it still evidences `migrantAge = 2`).

Likelihood given copy-type(s), per X-locus `l`:

- **Male** (one copy, type `tau`), allele `u`:
  ```
  L^X_l = psi_l(u | tau)
  ```
- **Female** (two copies, types `tau1, tau2`), alleles `(u,v)`:
  - both copies the same within-population type `p` (both `B`, or both `N`):
    ```
    L^X_l = H_p(u,v)          (apply inbreeding F_p)
    ```
  - mixed (`{B,N}`): outbred, sum over phase, no `F`:
    ```
    L^X_l = (1/2)[ f_{b,u} f_{s,v} + f_{b,v} f_{s,u} ]
    ```

Inbreeding is never applied to hemizygous males.

---

## 8. Prior over the latent state

Factor the state prior as BA3's existing ancestry prior times the new sex factor:

```
P(b, migrantAge, sigma | phi) = P0(migrantAge, b) * W(sigma | migrantAge, g_i, phi)
```

`P0(migrantAge, b)` is BA3's current ancestry prior (unchanged). The new factor
`W` reflects that the first-generation migrant is female with probability `phi`:

```
migrantAge = 0 : W = 1                                   (no migrant)

migrantAge = 1 : state consistent only if sigma = g_i
                 W = phi        if g_i = F
                 W = 1 - phi    if g_i = M

migrantAge = 2 : W(sigma = F) = phi
                 W(sigma = M) = 1 - phi
```

Only the first-generation migrant's sex is drawn from `phi`.

---

## 9. Marginal individual likelihood

```
L_i =   P0(0)      * A_i(0) * X_i(0)
      + sum_{b != s} P0(1, b) * W(g_i | 1, g_i, phi) * A_i(1, b) * X_i(1, b)
      + sum_{b != s} sum_{sigma in {F,M}}
            P0(2, b) * W(sigma | 2, g_i, phi) * A_i(2, b) * X_i(2, b, g_i, sigma)
```

where

```
A_i(.) = prod_{l in L_A} [ autosomal per-locus term, Section 7.1 ]
X_i(.) = prod_{l in L_X} [ X per-locus term,        Section 7.2 ]
```

The only new sum is over `sigma in {F, M}` for `migrantAge = 2`, and for females
the two `sigma` terms are identical (both give `{B,N}`), so they combine into a
single X evaluation weighted by `W(F) + W(M) = 1`. As in current BA3 the sampler
may instead **sample** the discrete state (Section 10).

---

## 10. MCMC updates

BA3 samples explicit per-individual migrant ancestry (`proposeMigrantAncDrop` /
`proposeMigrantAncAdd`). The extension adds the `sigma` label and one global
parameter update.

1. **Ancestry + sex moves.** For `migrantAge = 2`, extend the proposals to also
   draw `sigma in {F, M}`, scored by `W` and the X/autosomal likelihood. For
   `migrantAge = 1`, `sigma = g_i` is fixed. For females at `migrantAge = 2` the
   X likelihood does not depend on `sigma`, so `sigma` is updated by its prior
   (Gibbs from `phi`) there.

2. **Female fraction `phi` — conjugate Gibbs.** Each migrant-ancestry individual
   contributes one Bernoulli(`phi`) observation: the sex of its first-generation
   migrant (the focal's own sex at `migrantAge = 1`; the migrant parent's sex at
   `migrantAge = 2`). With a `Beta(a0, b0)` prior:
   ```
   phi | . ~ Beta(a0 + N_F, b0 + N_M)
   ```
   where `N_F` / `N_M` = number of individuals whose current first-generation
   migrant is female / male. Default `a0 = b0 = 1` (uniform).

3. **Baseline rates `M_{s<-b}`.** Unchanged from current BA3 — the sex split
   enters only through `W`, and the constraint `sum_b M_{s<-b} <= 1` is untouched.

4. **Savage-Dickey test for `phi = 1/2`.** Reuse the existing density-ratio
   machinery to report a Bayes factor for no sex bias.

Monitor acceptance rates for the ancestry/sex moves in the usual 20-60% window.

---

## 11. Identifiability

- `phi` is informed most directly by the **observed sexes of first-generation
  migrants** (`migrantAge = 1`); this works even with autosomes only.
- X-linked markers add power through **second-generation males**, whose whole-X
  ancestry (source vs. native) reads the migrant parent's sex outright.
- Second-generation **females** contribute no sex-of-migrant information from the
  X. Without X-linked loci, `phi` is identified solely from `migrantAge = 1`
  migrant sexes — usable but weaker; the X markers (via 2nd-gen males) sharpen it.

---

## 12. Generalization to more than two populations

The single global `phi` and all logic are unchanged: `b` sums over all sources
`!= s`, the native background frequencies are always `f_s`, and the copy types in
Section 7.2 apply with `B = f_b`, `N = f_s`. A later refinement could allow a
per-pair or directional `phi_{s<-b}`; the recommended starting point is one
global `phi`.

---

## 13. Implementation mapping (BA3 code)

- `include/BA3.h`
  - add `sex` to `struct indiv`; add a per-locus type annotation (autosomal / X)
    and hemizygous coding for male X;
  - add global `phi`; migration storage stays the baseline `M` matrix (no
    separate male/female rate matrices).
- `readMetadataFile` / `readVCFFile` : parse the sex column and the X-locus
  annotation; exclude PAR.
- `logLik` / `oneLocusLogLik` : branch on locus type and individual sex; add the
  X copy-type terms of Section 7.2 (autosomal path unchanged, still ages 0/1/2).
- `proposeMigrantAncDrop` / `proposeMigrantAncAdd` : carry the `sigma` label for
  `migrantAge = 2`.
- new `updatePhi` : the Beta-Bernoulli Gibbs step of Section 10.2.
- Savage-Dickey functions : add the `phi = 1/2` test.

---

## 14. Possible extension (not in scope)

Extending the ancestry depth to a third generation (a grandparent migrant) would
introduce a **recombinant X mosaic**: a third-generation individual whose bridge
parent is the mother inherits an X that is a per-marker mixture of source and
native (under the unlinked-loci assumption, each X marker source w.p. 1/2). That
would add power to estimate `phi` from third-generation individuals, at the cost
of a deeper ancestry state and a mosaic (`R`) copy-type in the X likelihood. It
is deliberately left out here to keep BA3's original two-generation depth.

---

## 15. Open items

- Default `Beta(a0, b0)` prior for `phi` (uniform vs weakly centered at 1/2).
- Sex inference from data (e.g., Y-linked or X-heterozygosity signal) as an
  alternative to a metadata sex column — out of scope for the first version.
