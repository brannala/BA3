#!/usr/bin/env python3
"""Exact BA3 posterior for a tiny data set, by enumeration (no MCMC).

The BA3 model:
  * migration rates: for each (receiving) population s, the row (m_sj, j != s) is
    uniform on {m_sj >= 0, sum_j m_sj <= 1/3};
  * ancestry of an individual sampled in s: age 0 (native) w.p. 1 - 3*sum_j m_sj,
    age 1 (migrant from j) w.p. m_sj, age 2 (child of a migrant from j) w.p. 2*m_sj;
  * allele frequencies ~ Dirichlet(1) per population and locus;
  * inbreeding F_q ~ Uniform(0,1) per population; an age-0/1 genotype drawn from
    population q is homozygous a/a w.p. F*p_a + (1-F)*p_a^2 and heterozygous a/b
    w.p. 2*(1-F)*p_a*p_b; an age-2 genotype has one (outbred) copy from the source
    and one from the native population, phase unknown;
  * missing genotypes contribute nothing.

Everything except ancestry is integrated out exactly:
  * m: with u = 3m the prior is a flat Dirichlet on the simplex, so the labelled
    ancestry assignment of population s has probability
        (P-1)! 3^-C 2^n2 n0! prod_j c_j! / (n0 + C + P - 1)!
    (c_j = migrants from j of age 1 or 2, C = sum_j c_j, n2 = age-2 count), and
    m | ancestry is Dirichlet, giving E[m_sj] and E[m_sj^2] in closed form;
  * allele frequencies: Dirichlet-multinomial on the gene copies, after summing
    over each homozygote's IBD indicator and each age-2 heterozygote's phase;
  * F: given everything else the integrand is a polynomial in each F_q of degree
    <= (number of genotypes), so Gauss-Legendre quadrature with enough nodes is
    exact (up to rounding).
The (2P-1)^N ancestry assignments are enumerated, so keep N small (~8 for P=2).

usage: tests/exact_posterior.py data.txt [--json out.json]
Output keys (populations by name, as in tests/check_exact.py): m[recv][src],
F[pop], freq pop/locus/allele, anc individual[srcpop,age].
"""
import argparse, itertools, json, math, sys
from collections import defaultdict
import numpy as np


def read_ba3(fn):
    """Native BA3 format: indiv popln locus allele1 allele2 ('0' = missing)."""
    geno, pop_of, loci, alleles = {}, {}, [], defaultdict(set)
    ind_order, pop_order = [], []
    for line in open(fn):
        t = line.split()
        if len(t) != 5:
            continue
        i, p, l, a, b = t
        if i not in pop_of:
            ind_order.append(i)
        if p not in pop_order:
            pop_order.append(p)
        pop_of[i] = p
        if l not in loci:
            loci.append(l)
        for x in (a, b):
            if x != "0":
                alleles[l].add(x)
        geno[(i, l)] = None if (a == "0" or b == "0") else (a, b)
    return ind_order, pop_of, pop_order, loci, alleles, geno


def ba3_pop_index(pop_order):
    # internal numbering only; results are reported by population name
    return {p: k for k, p in enumerate(pop_order)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data")
    ap.add_argument("--json", default="")
    ap.add_argument("--nodes", type=int, default=0, help="Gauss-Legendre nodes per F (default: exact)")
    a = ap.parse_args()

    inds, pop_of, pops, loci, allele_sets, geno = read_ba3(a.data)
    P, N = len(pops), len(inds)
    pidx = ba3_pop_index(pops)
    pname = {v: k for k, v in pidx.items()}
    allele_list = {l: sorted(allele_sets[l]) for l in loci}
    K = {l: len(allele_list[l]) for l in loci}
    aidx = {l: {x: k for k, x in enumerate(allele_list[l])} for l in loci}
    samp = [pidx[pop_of[i]] for i in inds]
    G = {(n, li): (None if geno[(i, l)] is None else (aidx[l][geno[(i, l)][0]], aidx[l][geno[(i, l)][1]]))
         for n, i in enumerate(inds) for li, l in enumerate(loci)}

    # Gauss-Legendre nodes on (0,1), exact for polynomial degree <= 2*nodes-1
    deg = N * len(loci)
    nodes = a.nodes or (deg // 2 + 2)
    x, w = np.polynomial.legendre.leggauss(nodes)
    Fn, Fw = (x + 1) / 2, w / 2
    grid = np.meshgrid(*([Fn] * P), indexing="ij")          # P-dim grid of F values
    gw = np.ones_like(grid[0])
    for q in range(P):
        gw = gw * np.meshgrid(*([Fw] * P), indexing="ij")[q]
    Fg = [g.ravel() for g in grid]
    gw = gw.ravel()
    logF = [np.log(f) for f in Fg]
    log1mF = [np.log1p(-f) for f in Fg]

    lf = [math.lgamma(k + 1) for k in range(4 * N + 20)]   # log k!

    # ancestry states for an individual sampled in s: (src, age)
    def states(s):
        return [(s, 0)] + [(j, age) for j in range(P) if j != s for age in (1, 2)]

    # per (individual, locus, state): list of latent options
    #   (copies [(pop, allele), ...], ibd_pop_or_-1, is_ibd, trial_pop_or_-1, log const)
    def options(n, li, st):
        g = G[(n, li)]
        if g is None:
            return [((), -1, 0, -1, 0.0)]
        a0, a1 = g
        src, age = st
        if age in (0, 1):
            q = src                                   # copies come from src (= s for age 0)
            if a0 == a1:
                return [(((q, a0),), q, 1, q, 0.0),                  # IBD: one copy, factor F
                        (((q, a0), (q, a0)), q, 0, q, 0.0)]          # outbred: two copies, (1-F)
            return [(((q, a0), (q, a1)), q, 0, q, math.log(2.0))]    # het: 2(1-F) p_a p_b
        s = samp[n]
        if a0 == a1:
            return [(((src, a0), (s, a0)), -1, 0, -1, 0.0)]
        return [(((src, a0), (s, a1)), -1, 0, -1, 0.0),              # phase 1
                (((src, a1), (s, a0)), -1, 0, -1, 0.0)]              # phase 2

    St = [states(samp[n]) for n in range(N)]
    OPT = {(n, li, k): options(n, li, st) for n in range(N) for li in range(len(loci))
           for k, st in enumerate(St[n])}

    # accumulators (unnormalised, relative to a running log offset)
    logZ = -math.inf
    acc = defaultdict(float)          # key -> sum of weight * E[.]
    offset = None

    def add(key, logw, val):
        acc[key] += math.exp(logw - offset) * val

    configs = list(itertools.product(*[range(len(St[n])) for n in range(N)]))
    # first pass would be needed for a stable offset; use the first config's weight
    for ci, cfg in enumerate(configs):
        anc = [St[n][cfg[n]] for n in range(N)]
        # ---- prior of the labelled ancestry assignment, m integrated out
        logPA = 0.0
        Em, Em2 = {}, {}
        for s in range(P):
            members = [anc[n] for n in range(N) if samp[n] == s]
            n0 = sum(1 for src, age in members if age == 0)
            c = {j: sum(1 for src, age in members if src == j and age > 0) for j in range(P) if j != s}
            n2 = sum(1 for src, age in members if age == 2)
            C = sum(c.values())
            logPA += (lf[P - 1] - C * math.log(3.0) + n2 * math.log(2.0) + lf[n0]
                      + sum(lf[cj] for cj in c.values()) - lf[n0 + C + P - 1])
            a0 = n0 + C + P                                  # Dirichlet total
            for j, cj in c.items():
                Em[(s, j)] = (cj + 1) / a0 / 3.0
                Em2[(s, j)] = (cj + 1) * (cj + 2) / (a0 * (a0 + 1)) / 9.0
        # ---- genotype likelihood: per locus, sum over latent options, as a function of F
        # per locus: dict exponent-key -> scalar sum, plus freq moment sums
        loc_W, loc_M = [], []
        for li, l in enumerate(loci):
            Kl = K[l]
            W = defaultdict(float)                           # (kF tuple, trials tuple) -> sum
            M = defaultdict(lambda: np.zeros((P, Kl, 2)))    # same key -> sum of weight*(E p, E p^2)
            opts = [OPT[(n, li, cfg[n])] for n in range(N)]
            for combo in itertools.product(*opts):
                cnt = [[0] * Kl for _ in range(P)]
                kF = [0] * P
                tr = [0] * P
                lc = 0.0
                for copies, ibdp, isibd, trp, c0 in combo:
                    for q, al in copies:
                        cnt[q][al] += 1
                    if trp >= 0:
                        tr[trp] += 1
                        kF[trp] += isibd
                    lc += c0
                # Dirichlet(1)-multinomial of the ordered copies, per population
                ldm = 0.0
                for q in range(P):
                    tot = sum(cnt[q])
                    ldm += lf[Kl - 1] + sum(lf[v] for v in cnt[q]) - lf[tot + Kl - 1]
                val = math.exp(lc + ldm)
                key = (tuple(kF), tuple(tr))
                W[key] += val
                mom = M[key]
                for q in range(P):
                    tot = sum(cnt[q])
                    for al in range(Kl):
                        mom[q, al, 0] += val * (cnt[q][al] + 1) / (tot + Kl)
                        mom[q, al, 1] += val * (cnt[q][al] + 1) * (cnt[q][al] + 2) / ((tot + Kl) * (tot + Kl + 1))
            # evaluate each locus factor on the F grid: sum_key W * prod_q F^k (1-F)^(t-k)
            fac = np.zeros_like(gw)
            mfac = np.zeros((P, Kl, 2) + gw.shape)
            for key, v in W.items():
                kF, tr = key
                e = np.zeros_like(gw)
                for q in range(P):
                    e = e + kF[q] * logF[q] + (tr[q] - kF[q]) * log1mF[q]
                ev = np.exp(e)
                fac += v * ev
                mfac += M[key][..., None] * ev
            loc_W.append(fac)
            loc_M.append(mfac)
        prodW = np.prod(loc_W, axis=0)                       # P(G | anc, F) on the grid
        like = float(np.dot(gw, prodW))
        if like <= 0:
            continue
        logw = logPA + math.log(like)
        if offset is None:
            offset = logw
        if logw - offset > 600:                              # rescale accumulators
            sc = math.exp(offset - logw)
            for k in acc:
                acc[k] *= sc
            offset = logw
        wt = math.exp(logw - offset)
        acc["Z"] += wt
        for (s, j), v in Em.items():
            acc[("m", s, j, 1)] += wt * v
            acc[("m", s, j, 2)] += wt * Em2[(s, j)]
        for q in range(P):
            acc[("F", q, 1)] += wt * float(np.dot(gw, prodW * Fg[q])) / like
            acc[("F", q, 2)] += wt * float(np.dot(gw, prodW * Fg[q] ** 2)) / like
        for n in range(N):
            acc[("anc", n, anc[n])] += wt
        for li in range(len(loci)):
            others = prodW / np.where(loc_W[li] == 0, 1, loc_W[li])
            for q in range(P):
                for al in range(K[loci[li]]):
                    for mom in (0, 1):
                        acc[("freq", q, li, al, mom + 1)] += wt * float(
                            np.dot(gw, loc_M[li][q, al, mom] * others)) / like

    Z = acc["Z"]
    out = {}
    for s in range(P):
        for j in range(P):
            if s == j:
                continue
            m1, m2 = acc[("m", s, j, 1)] / Z, acc[("m", s, j, 2)] / Z
            out[f"m[{pname[s]}][{pname[j]}]"] = (m1, math.sqrt(max(m2 - m1 * m1, 0)))
    for q in range(P):
        f1, f2 = acc[("F", q, 1)] / Z, acc[("F", q, 2)] / Z
        out[f"F[{pname[q]}]"] = (f1, math.sqrt(max(f2 - f1 * f1, 0)))
    for li, l in enumerate(loci):
        for q in range(P):
            for al, name in enumerate(allele_list[l]):
                f1 = acc[("freq", q, li, al, 1)] / Z
                f2 = acc[("freq", q, li, al, 2)] / Z
                out[f"freq {pname[q]}/{l}/{name}"] = (f1, math.sqrt(max(f2 - f1 * f1, 0)))
    for n, i in enumerate(inds):
        for st in St[n]:
            out[f"anc {i}[{pname[st[0]]},{st[1]}]"] = (acc[("anc", n, st)] / Z, None)
    if a.json:
        json.dump({"post": out}, open(a.json, "w"), indent=1)
    for k in sorted(out):
        if not k.startswith("anc") and not k.startswith("freq"):
            print(f"{k:10s} mean {out[k][0]:.5f}  SD {out[k][1]:.5f}")
    print(f"({len(configs)} ancestry configurations, {nodes} quadrature nodes per F)")


if __name__ == "__main__":
    main()
