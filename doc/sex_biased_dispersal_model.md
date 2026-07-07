# Sex-Biased Dispersal Model for BA3

**Status:** design draft
**Scope:** extends BA3's migrant-ancestry model to infer a global dispersal
sex-bias from autosomal + X-linked markers with per-individual sex metadata,
keeping BA3's existing two-generation ancestry depth (`migrantAge` 0, 1, 2).

---

## 1. Goal

BA3 currently estimates a single directional migration rate `M_{s<-b}` (fraction
of population `s` that is migrant-derived from source `b` per generation). This
extension infers the **sex bias of dispersal** from autosomal + X-linked markers
with per-individual sex, distinguishing three global female fractions — movement
(`phiMove`), effective gene flow (`phiBreed`), and the resident baseline (`rho`)
— while retaining all standard BA3 assumptions (HWE within populations, constant
allele frequencies, unlinked loci, inbreeding coefficient `F` per population) and
BA3's existing ancestry depth.

The signal comes from **X-linked inheritance being sex-asymmetric**: a male
transmits his single X only to daughters, and a second-generation individual's X
therefore records whether its migrant parent was the mother or the father. This
lets the data inform the sex composition of migrants — and, because a
second-generation individual exists only if its migrant parent *bred*, it
separates who *moves* from who actually contributes *genes* (Section 3).

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

## 3. Parametrization (female fractions: movement, gene flow, residents)

The sex bias of dispersal is described by **global female fractions** rather than
a fixed 50:50 assumption. Three quantities are inferred:

- `phiMove in (0,1)`  = female fraction of **age-1 migrants** (individuals that
  **moved**; a mover may or may not breed) — *movement dispersal*.
- `phiBreed in (0,1)` = female fraction of **age-2 migrant parents** (migrants
  that moved **and bred**, since an age-2 individual only exists because its
  migrant parent reproduced) — *effective gene flow*.
- `rho in (0,1)`      = female fraction of **residents** (non-dispersers: the own
  sex of age-0 and age-2 individuals) — the *population sex ratio* baseline.

**Movement vs gene flow.** `phiMove` and `phiBreed` are distinct because
migration is not gene flow: a first-generation migrant that fails to breed (e.g.
a male excluded by reproductive competition / high male reproductive variance)
contributes movement but no genes. `phiBreed` is the population-genetically
meaningful quantity, and — because a male's single X comes from his mother — the
**X-linked data specifically informs `phiBreed`** (via the age-2 male
migrant-parent sex). `phiMove` is informed by the observed sexes of age-1
migrants.

**Bias is measured against `rho`, not 1/2.** Comparing a migrant fraction to a
fixed 1/2 conflates dispersal sex bias with the population/sampling sex ratio.
The dispersal sex bias is therefore `phi vs rho`: `phi > rho` female-biased,
`phi < rho` male-biased, `phi = rho` no bias. Inferring `rho` from non-migrants
controls for a skewed population or sex-biased sampling (as long as the skew is
the same for residents and migrants).

The baseline migration rate `M_{s<-b}` retains its usual meaning (sex-summed);
`phiMove`/`phiBreed` split the migrant sexes at each generation. Testing
`phiMove = rho`, `phiBreed = rho`, and `phiMove = phiBreed` (does breeding
success differ from movement?) are Savage-Dickey density-ratio tests, analogous
to BA3's existing `M = 0` tests.

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

1. **Ancestry + sex moves.** For `migrantAge = 2`, the proposal also draws
   `sigma in {F, M}` from its prior `Bernoulli(phiBreed)` (so the sigma prior and
   proposal densities cancel in the MH ratio), scored by the X/autosomal
   likelihood. For `migrantAge = 1`, `sigma = g_i` is fixed. The residual per-
   individual sex term added to the MH ratio is `log(phiMove)` (age-1 own sex) or
   `log(rho)` (age-0 / age-2 own sex).

2. **Female fractions — conjugate Gibbs.** With a `Beta(a0, b0)` prior
   (`a0 = b0 = 1` uniform), each is updated from Bernoulli counts of the relevant
   sexes:
   ```
   phiMove  | . ~ Beta(a0 + age1_F,   b0 + age1_M)     # age-1 migrant own sexes
   phiBreed | . ~ Beta(a0 + age2sig_F, b0 + age2sig_M)  # age-2 migrant-parent sexes
   rho      | . ~ Beta(a0 + res_F,     b0 + res_M)      # age-0 + age-2 own sexes
   ```
   The age-2 `sigma` Gibbs step uses `phiBreed` (not `phiMove`) as its prior.

3. **Baseline rates `M_{s<-b}`.** Unchanged from current BA3.

4. **Savage-Dickey tests (Rao-Blackwellized).** Because each fraction has a
   conjugate Beta full-conditional, the posterior density needed for each test is
   the average of the conjugate density over MCMC samples (more accurate than a
   KDE). The density of a difference at 0 is the overlap integral of the two
   Beta full-conditionals. Reported nulls: `phiMove = rho` (sex-biased movement),
   `phiBreed = rho` (sex-biased gene flow), `phiMove = phiBreed` (movement differs
   from gene flow).

Monitor acceptance rates for the ancestry/sex moves in the usual 20-60% window.

---

## 11. Identifiability

- `phiMove` is informed by the **observed sexes of first-generation migrants**
  (`migrantAge = 1`); it needs no X data.
- `phiBreed` is informed by **second-generation males**, whose whole-X ancestry
  (source vs. native) reads the migrant parent's sex outright — so `phiBreed`
  effectively requires X-linked loci. Second-generation **females** contribute no
  sex-of-migrant information from the X.
- `rho` is informed by the observed sexes of residents (age-0 and age-2).
- **Caveat (X/autosome discordance).** Age classification (age-0/1/2) uses the
  whole genome including the X. If X and autosomal ancestry are discordant for
  reasons *other* than recent sex-biased gene flow (e.g. deep sex-specific
  demography, lower X `Ne`), that discordance can be misread as `phiBreed` bias
  and can even shift which individuals are called age-1 vs age-2 (contaminating
  `phiMove`). Comparing the full model to an **autosome-only** run isolates the
  movement signal and diagnoses such discordance.

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
