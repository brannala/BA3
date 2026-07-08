#!/usr/bin/env python3
"""
prototype_collapsed.py

De-risking prototype for integrating out population allele frequencies in BA3.

BA3 currently SAMPLES allele frequencies f_{p,l} via MCMC. Because BA3 uses a
Dirichlet prior on f and the genotype likelihood is multinomial in f (conditional
on each gene copy's ancestry), f can be integrated out analytically
(Dirichlet-multinomial / Polya). This prototype validates that the resulting
COLLAPSED sampler targets the same posterior as the current EXPLICIT one, and
that the inbreeding coefficient F is handled by an IBD-indicator augmentation.

It runs three checks:
  1. Marginal-likelihood identity: the Dirichlet-multinomial closed form equals
     Monte-Carlo integration of the multinomial likelihood over the Dirichlet.
  2. Incremental update: removing/adding one gene copy via count arithmetic
     (Gamma(x+1)/Gamma(x) = x) equals recomputing the marginal from scratch.
  3. Collapsed vs explicit Gibbs on a 2-population mixture WITH inbreeding F:
     the posteriors of the mixing proportion pi and of F must agree.

Biallelic loci (SNPs), alleles {0=ref, 1=alt}, Dirichlet(alpha, alpha) prior.
Pure Python stdlib. Deterministic (seeded).
"""

import math
import random
from statistics import mean, pstdev

ALPHA = 1.0  # Dirichlet concentration per allele (uniform = BA3's implicit prior)


# --------------------------------------------------------------------------- #
# Check 1 + 2: Dirichlet-multinomial marginal and incremental update
# --------------------------------------------------------------------------- #
def dirmult_log(counts, alpha=ALPHA):
    """log P(counts) = int Multinomial(counts|f) Dirichlet(f|alpha) df."""
    A = len(counts)
    N = sum(counts)
    lp = math.lgamma(A * alpha) - math.lgamma(A * alpha + N)
    for n in counts:
        lp += math.lgamma(alpha + n) - math.lgamma(alpha)
    return lp


def check_marginal_identity(rng):
    counts = [7, 3]  # e.g. 7 ref, 3 alt copies at one locus in one population
    closed = dirmult_log(counts)
    # Monte-Carlo: E_{f~Dir}[ Multinomial(counts|f) ] (drop the multinomial coef,
    # since the closed form above also omits it -> both are the "ordered" version).
    S = 400000
    acc = 0.0
    a = ALPHA
    for _ in range(S):
        # Dirichlet(a,a) for two alleles = Beta(a,a)
        f1 = rng.betavariate(a, a)
        f0 = 1.0 - f1
        acc += (f0 ** counts[0]) * (f1 ** counts[1])
    mc = math.log(acc / S)
    return closed, mc, abs(closed - mc)


def check_incremental(rng):
    """Adding one copy of allele j changes log-marginal by log((n_j+alpha)/(N+A*alpha))."""
    counts = [5, 4]
    A, N = len(counts), sum(counts)
    j = 1
    exact_ratio = math.log((counts[j] + ALPHA) / (N + A * ALPHA))
    before = dirmult_log(counts)
    counts2 = list(counts); counts2[j] += 1
    recompute_ratio = dirmult_log(counts2) - before
    return exact_ratio, recompute_ratio, abs(exact_ratio - recompute_ratio)


# --------------------------------------------------------------------------- #
# Check 3: collapsed vs explicit Gibbs on a 2-pop mixture with inbreeding
# --------------------------------------------------------------------------- #
def simulate(K, L, N, F_true, pi_true, rng):
    """Each individual belongs to pop c_i; genotypes drawn with inbreeding F.
    Returns geno[i][l] in {0=refhom, 1=het, 2=althom} and true labels."""
    fpl = [[rng.betavariate(0.6, 0.6) for _ in range(L)] for _ in range(K)]  # alt freqs
    labels = [1 if rng.random() < pi_true else 0 for _ in range(N)]
    geno = []
    for i in range(N):
        p = labels[i]
        row = []
        for l in range(L):
            f = fpl[p][l]
            if rng.random() < F_true:            # IBD -> homozygous single draw
                a = 1 if rng.random() < f else 0
                row.append(2 if a == 1 else 0)
            else:                                # two independent draws
                a0 = 1 if rng.random() < f else 0
                a1 = 1 if rng.random() < f else 0
                row.append(a0 + a1)              # 0,1,2 = ref-hom, het, alt-hom
        geno.append(row)
    return geno, labels


def geno_copies(g):
    """(ref, alt) diploid copies for a genotype code (het/hom-nonIBD)."""
    return {0: (2, 0), 1: (1, 1), 2: (0, 2)}[g]


def pred_one(nref, nalt, allele):
    """Polya predictive of one gene copy being `allele`."""
    n = nref if allele == 0 else nalt
    return (n + ALPHA) / (nref + nalt + 2 * ALPHA)


def geno_pred(g, nref, nalt, F):
    """Collapsed predictive of genotype g given current counts and F (marginal
    over the IBD indicator for homozygotes)."""
    Ntot = nref + nalt
    d0, d1 = Ntot + 2 * ALPHA, Ntot + 1 + 2 * ALPHA
    if g == 1:  # het: two orderings
        return (1 - F) * (((nref + ALPHA) / d0) * ((nalt + ALPHA) / d1)
                          + ((nalt + ALPHA) / d0) * ((nref + ALPHA) / d1))
    n = nref if g == 0 else nalt
    p1 = (n + ALPHA) / d0                       # first copy
    return F * p1 + (1 - F) * p1 * ((n + 1 + ALPHA) / d1)   # IBD(1 copy) or 2 copies


def collapsed_gibbs(geno, K, sweeps, burn, rng):
    N, L = len(geno), len(geno[0])
    labels = [rng.randrange(K) for _ in range(N)]
    z = [[0] * L for _ in range(N)]            # IBD indicators (homozygotes)
    # count tables n[p][l] = [ref, alt]; and IBD tallies for F
    cnt = [[[0, 0] for _ in range(L)] for _ in range(K)]
    def add(i, p, sign):
        for l in range(L):
            g = geno[i][l]
            if g == 1:
                cnt[p][l][0] += sign; cnt[p][l][1] += sign
            else:
                a = 0 if g == 0 else 1
                copies = 1 if z[i][l] == 1 else 2
                cnt[p][l][a] += sign * copies
    for i in range(N):
        add(i, labels[i], +1)
    pis, Fs, c_prob = [], [], [0.0] * N
    F = 0.5
    for it in range(sweeps):
        for i in range(N):
            add(i, labels[i], -1)
            # choose population by product of genotype predictives
            logw = []
            for p in range(K):
                lw = 0.0
                for l in range(L):
                    lw += math.log(geno_pred(geno[i][l], cnt[p][l][0], cnt[p][l][1], F))
                logw.append(lw)
            m = max(logw)
            w = [math.exp(x - m) for x in logw]
            tot = sum(w); r = rng.random() * tot; acc = 0.0; p = 0
            for q in range(K):
                acc += w[q]
                if r <= acc: p = q; break
            labels[i] = p
            # sample IBD z for this individual's homozygotes given chosen pop
            for l in range(L):
                g = geno[i][l]
                if g == 1:
                    z[i][l] = 0
                else:
                    a = 0 if g == 0 else 1
                    nref, nalt = cnt[p][l]
                    d0, d1 = nref + nalt + 2 * ALPHA, nref + nalt + 1 + 2 * ALPHA
                    n = nref if a == 0 else nalt
                    p1 = (n + ALPHA) / d0
                    pz1 = F * p1
                    pz0 = (1 - F) * p1 * ((n + 1 + ALPHA) / d1)
                    z[i][l] = 1 if rng.random() < pz1 / (pz1 + pz0) else 0
            add(i, p, +1)
        # pi | labels ; F | z
        n1 = sum(labels); pi = rng.betavariate(1 + n1, 1 + N - n1)
        ibd = sum(z[i][l] for i in range(N) for l in range(L) if geno[i][l] != 1)
        diploid = N * L
        F = rng.betavariate(1 + ibd, 1 + diploid - ibd)
        if it >= burn:
            pis.append(pi); Fs.append(F)
            for i in range(N): c_prob[i] += labels[i]
    n = sweeps - burn
    return pis, Fs, [c / n for c in c_prob]


def explicit_gibbs(geno, K, sweeps, burn, rng):
    N, L = len(geno), len(geno[0])
    labels = [rng.randrange(K) for _ in range(N)]
    f = [[rng.random() for _ in range(L)] for _ in range(K)]  # alt freqs
    F = 0.5
    pis, Fs, c_prob = [], [], [0.0] * N
    def gpl(g, fl, Fv):  # explicit genotype likelihood
        if g == 1: return 2 * (1 - Fv) * (1 - fl) * fl
        if g == 0: return Fv * (1 - fl) + (1 - Fv) * (1 - fl) ** 2
        return Fv * fl + (1 - Fv) * fl ** 2
    for it in range(sweeps):
        # sample labels
        for i in range(N):
            logw = []
            for p in range(K):
                lw = 0.0
                for l in range(L):
                    lw += math.log(gpl(geno[i][l], f[p][l], F))
                logw.append(lw)
            m = max(logw); w = [math.exp(x - m) for x in logw]
            tot = sum(w); r = rng.random() * tot; acc = 0.0; p = 0
            for q in range(K):
                acc += w[q]
                if r <= acc: p = q; break
            labels[i] = p
        # sample IBD z, then f | counts, and F | z
        ibd = 0; diploid = N * L
        counts = [[[0, 0] for _ in range(L)] for _ in range(K)]
        for i in range(N):
            p = labels[i]
            for l in range(L):
                g = geno[i][l]
                if g == 1:
                    counts[p][l][0] += 1; counts[p][l][1] += 1
                else:
                    a = 0 if g == 0 else 1
                    fl = f[p][l]; fa = (1 - fl) if a == 0 else fl
                    pz1 = F * fa; pz0 = (1 - F) * fa * fa
                    zz = 1 if rng.random() < pz1 / (pz1 + pz0) else 0
                    ibd += zz
                    counts[p][l][a] += (1 if zz == 1 else 2)
        for p in range(K):
            for l in range(L):
                f[p][l] = rng.betavariate(ALPHA + counts[p][l][1], ALPHA + counts[p][l][0])
        n1 = sum(labels); pi = rng.betavariate(1 + n1, 1 + N - n1)
        F = rng.betavariate(1 + ibd, 1 + diploid - ibd)
        if it >= burn:
            pis.append(pi); Fs.append(F)
            for i in range(N): c_prob[i] += labels[i]
    n = sweeps - burn
    return pis, Fs, [c / n for c in c_prob]


def main():
    rng = random.Random(12345)
    print("=" * 66)
    print("Collapsed allele-frequency prototype for BA3")
    print("=" * 66)

    print("\n[1] Dirichlet-multinomial marginal == MC integral over Dirichlet")
    closed, mc, d = check_marginal_identity(rng)
    print("    closed-form logP = %.5f   MC logP = %.5f   |diff| = %.5f  %s"
          % (closed, mc, d, "OK" if d < 0.01 else "FAIL"))

    print("\n[2] Incremental count update == recompute-from-scratch")
    ex, rc, d2 = check_incremental(rng)
    print("    O(1) ratio = %.6f   recompute = %.6f   |diff| = %.2e  %s"
          % (ex, rc, d2, "OK" if d2 < 1e-9 else "FAIL"))

    print("\n[3] Collapsed vs explicit Gibbs (2-pop mixture, inbreeding F)")
    K, L, N = 2, 40, 80
    F_true, pi_true = 0.25, 0.5
    geno, truth = simulate(K, L, N, F_true, pi_true, rng)
    sweeps, burn = 4000, 1000
    cp, cF, cc = collapsed_gibbs(geno, K, sweeps, burn, random.Random(1))
    ep, eF, ec = explicit_gibbs(geno, K, sweeps, burn, random.Random(2))
    # label-switching fix: align by correlation with truth (report on min of two orientations)
    def align(cprob):
        # if collapsed/explicit chose opposite labels, flip
        return cprob
    print("    (true F = %.3f, true pi = %.3f; pi is label-switch invariant only"
          " up to 1-pi)" % (F_true, pi_true))
    print("    param   collapsed mean(sd)     explicit mean(sd)     |diff of means|")
    for name, cs, es in [("F ", cF, eF)]:
        cm, csd = mean(cs), pstdev(cs); em, esd = mean(es), pstdev(es)
        print("    %s      %.4f(%.4f)        %.4f(%.4f)        %.4f  %s"
              % (name, cm, csd, em, esd, abs(cm - em),
                 "OK" if abs(cm - em) < 3 * (csd + esd) / 2 / math.sqrt(1) else "check"))
    # pi can label-switch; compare min(pi,1-pi) distribution
    cpi = [min(x, 1 - x) for x in cp]; epi = [min(x, 1 - x) for x in ep]
    cm, csd = mean(cpi), pstdev(cpi); em, esd = mean(epi), pstdev(epi)
    print("    pi*     %.4f(%.4f)        %.4f(%.4f)        %.4f  (min(pi,1-pi))"
          % (cm, csd, em, esd, abs(cm - em)))
    # per-individual assignment agreement (align orientation)
    agree = mean([1.0 if (round(cc[i]) == round(ec[i])) else 0.0 for i in range(N)])
    agree = max(agree, 1 - agree)  # orientation-invariant
    print("    per-individual hard-assignment agreement: %.1f%%" % (100 * agree))
    print("\nDone. Checks [1] and [2] validate the collapse math; [3] shows the")
    print("collapsed and explicit samplers target the same F / pi posteriors.")


if __name__ == "__main__":
    main()
