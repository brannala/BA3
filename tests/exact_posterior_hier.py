#!/usr/bin/env python3
"""Exact BA3 posterior under the hierarchical (F-model) allele-frequency prior.

Same model as tests/exact_posterior.py except for the allele frequencies:
  * ancestral frequencies pi_l ~ Dirichlet(1) per locus, shared by populations;
  * p_ql | pi_l ~ Dirichlet(c_q pi_l), c_q = (1 - theta_q) / theta_q;
  * theta_q ~ Uniform(0, 1) per population.
p is integrated out analytically (Dirichlet-multinomial with alpha = c_q pi_l);
pi_l by quadrature on the simplex (stick-breaking, Gauss-Legendre per stick);
theta and F by Gauss-Legendre quadrature on (0,1)^P. F is exact as before; pi
and theta are smooth but not polynomial, so check convergence with --nodes-pi
and --nodes-theta. The (2P-1)^N ancestry assignments are enumerated.

usage: tests/exact_posterior_hier.py data.txt [--json out.json]
Output keys as in exact_posterior.py, plus theta[pop] and pi locus/allele.
"""
import argparse, itertools, json, math, os, sys
from collections import defaultdict
import numpy as np
gammaln = np.vectorize(math.lgamma, otypes=[float])

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_posterior import read_ba3, ba3_pop_index  # noqa: E402


def gl01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return (x + 1) / 2, w / 2


def simplex_rule(K, n):
    """Nodes (M x K) and weights (M) for the uniform Dirichlet(1) density on the
    K-simplex, by stick-breaking: pi_1 = u_1, pi_k = u_k prod_{j<k}(1-u_j).
    Weights include the density (K-1)! and the Jacobian, so they sum to 1."""
    if K == 1:
        return np.ones((1, 1)), np.ones(1)
    u, w = gl01(n)
    pts, wts = [], []
    for idx in itertools.product(range(n), repeat=K - 1):
        # Jacobian of the stick-breaking map is prod_k rest_k, where rest_k is
        # the stick left before break k.
        rest, p, wt = 1.0, [], math.factorial(K - 1)
        for i in idx:
            p.append(rest * u[i])
            wt *= w[i] * rest
            rest *= 1 - u[i]
        p.append(rest)
        pts.append(p)
        wts.append(wt)
    return np.array(pts), np.array(wts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data")
    ap.add_argument("--json", default="")
    ap.add_argument("--nodes", type=int, default=0, help="Gauss-Legendre nodes per F (default: exact)")
    ap.add_argument("--nodes-theta", type=int, default=16)
    ap.add_argument("--nodes-pi", type=int, default=16, help="nodes per stick of the pi simplex")
    a = ap.parse_args()

    inds, pop_of, pops, loci, allele_sets, geno = read_ba3(a.data)
    P, N, L = len(pops), len(inds), len(loci)
    pidx = ba3_pop_index(pops)
    pname = {v: k for k, v in pidx.items()}
    allele_list = {l: sorted(allele_sets[l]) for l in loci}
    K = [len(allele_list[l]) for l in loci]
    aidx = {l: {x: k for k, x in enumerate(allele_list[l])} for l in loci}
    samp = [pidx[pop_of[i]] for i in inds]
    G = {(n, li): (None if geno[(i, l)] is None else (aidx[l][geno[(i, l)][0]], aidx[l][geno[(i, l)][1]]))
         for n, i in enumerate(inds) for li, l in enumerate(loci)}

    # ---- quadrature grids. Joint grid over (F_0..F_{P-1}, theta_0..theta_{P-1}),
    # flattened; F and theta axes are kept separate so terms can be outer products.
    deg = N * L
    nF = a.nodes or (deg // 2 + 2)
    Fn, Fw = gl01(nF)
    Tn, Tw = gl01(a.nodes_theta)
    Fg = [g.ravel() for g in np.meshgrid(*([Fn] * P), indexing="ij")]
    Fwt = np.prod([g.ravel() for g in np.meshgrid(*([Fw] * P), indexing="ij")], axis=0)
    Tg = [g.ravel() for g in np.meshgrid(*([Tn] * P), indexing="ij")]
    Twt = np.prod([g.ravel() for g in np.meshgrid(*([Tw] * P), indexing="ij")], axis=0)
    logF = [np.log(f) for f in Fg]
    log1mF = [np.log1p(-f) for f in Fg]
    cT = (1 - Tn) / Tn                                   # c on the 1-D theta nodes
    PI = [simplex_rule(K[li], a.nodes_pi) for li in range(L)]

    lf = [math.lgamma(k + 1) for k in range(4 * N + 20)]

    def states(s):
        return [(s, 0)] + [(j, age) for j in range(P) if j != s for age in (1, 2)]

    def options(n, li, st):
        g = G[(n, li)]
        if g is None:
            return [((), 0, -1, 0.0)]
        a0, a1 = g
        src, age = st
        if age in (0, 1):
            q = src
            if a0 == a1:
                return [(((q, a0),), 1, q, 0.0), (((q, a0), (q, a0)), 0, q, 0.0)]
            return [(((q, a0), (q, a1)), 0, q, math.log(2.0))]
        s = samp[n]
        if a0 == a1:
            return [(((src, a0), (s, a0)), 0, -1, 0.0)]
        return [(((src, a0), (s, a1)), 0, -1, 0.0), (((src, a1), (s, a0)), 0, -1, 0.0)]

    St = [states(samp[n]) for n in range(N)]
    OPT = {(n, li, k): options(n, li, st) for n in range(N) for li in range(L)
           for k, st in enumerate(St[n])}

    # ---- memoised building blocks
    fpoly_memo = {}

    def fpoly(kF, tr):
        """prod_q F_q^k (1-F_q)^(t-k) on the F grid."""
        key = (kF, tr)
        if key not in fpoly_memo:
            e = np.zeros_like(Fwt)
            for q in range(P):
                e = e + kF[q] * logF[q] + (tr[q] - kF[q]) * log1mF[q]
            fpoly_memo[key] = np.exp(e)
        return fpoly_memo[key]

    dm_memo = {}

    def dm_terms(li, cnt):
        """For count table cnt (P x K), return arrays over the theta grid:
        like[t]       = int prod_q DM(cnt_q | c_q pi) dpi
        pimom[a,t]    = int pi_a * prod_q DM dpi
        pmom[q,a,m,t] = int E[p_qa^m | ...] * prod_q DM dpi   (m = 1, 2)
        DM here is the ordered-sequence Dirichlet-multinomial probability."""
        key = (li, cnt)
        if key in dm_memo:
            return dm_memo[key]
        pts, wts = PI[li]                               # (M x Kl), (M)
        Kl = K[li]
        arr = np.array(cnt, dtype=float)               # P x Kl
        # per population: log DM over (M pi nodes) x (nT theta nodes)
        per_q, alpha_q = [], []
        for q in range(P):
            al = cT[None, :, None] * pts[:, None, :]   # M x nT x Kl
            n = arr[q][None, None, :]
            Nq = arr[q].sum()
            ld = (gammaln(cT)[None, :] - gammaln(cT + Nq)[None, :]
                  + (gammaln(al + n) - gammaln(al)).sum(axis=2))
            per_q.append(ld)
            alpha_q.append(al)
        # joint over theta grid: sum of per-q logs broadcast to M x nT x ... x nT
        shape = [len(wts)] + [len(Tn)] * P
        tot = np.zeros(shape)
        for q in range(P):
            sh = [len(wts)] + [1] * P
            sh[1 + q] = len(Tn)
            tot = tot + per_q[q].reshape(sh)
        lik = np.exp(tot)                              # M x nT^P
        wl = lik * wts.reshape([-1] + [1] * P)
        like = wl.sum(axis=0).ravel()
        pimom = np.stack([(wl * pts[:, k].reshape([-1] + [1] * P)).sum(axis=0).ravel()
                          for k in range(Kl)])
        pmom = np.zeros((P, Kl, 2, len(Twt)))
        for q in range(P):
            Nq = arr[q].sum()
            sh = [len(wts)] + [1] * P
            sh[1 + q] = len(Tn)
            for k in range(Kl):
                num = alpha_q[q][:, :, k] + arr[q][k]              # M x nT
                den = cT[None, :] + Nq
                e1 = (num / den).reshape(sh)
                e2 = (num * (num + 1) / (den * (den + 1))).reshape(sh)
                pmom[q, k, 0] = (wl * e1).sum(axis=0).ravel()
                pmom[q, k, 1] = (wl * e2).sum(axis=0).ravel()
        dm_memo[key] = (like, pimom, pmom)
        return dm_memo[key]

    acc = defaultdict(float)
    offset = None
    configs = list(itertools.product(*[range(len(St[n])) for n in range(N)]))
    for cfg in configs:
        anc = [St[n][cfg[n]] for n in range(N)]
        logPA = 0.0
        Em, Em2, Z0 = {}, {}, {}
        for s in range(P):
            members = [anc[n] for n in range(N) if samp[n] == s]
            n0 = sum(1 for src, age in members if age == 0)
            c = {j: sum(1 for src, age in members if src == j and age > 0) for j in range(P) if j != s}
            n2 = sum(1 for src, age in members if age == 2)
            C = sum(c.values())
            logPA += (lf[P - 1] - C * math.log(3.0) + n2 * math.log(2.0) + lf[n0]
                      + sum(lf[cj] for cj in c.values()) - lf[n0 + C + P - 1])
            a0 = n0 + C + P
            for j, cj in c.items():
                Em[(s, j)] = (cj + 1) / a0 / 3.0
                Em2[(s, j)] = (cj + 1) * (cj + 2) / (a0 * (a0 + 1)) / 9.0
                Z0[(s, j)] = cj == 0

        # per locus: factor over (F grid) x (theta grid), and its moment versions
        loc = []
        for li in range(L):
            Kl = K[li]
            groups = defaultdict(lambda: defaultdict(float))   # cnt -> (kF,tr) -> weight
            opts = [OPT[(n, li, cfg[n])] for n in range(N)]
            for combo in itertools.product(*opts):
                cnt = [[0] * Kl for _ in range(P)]
                kF, tr = [0] * P, [0] * P
                lc = 0.0
                for copies, isibd, trp, c0 in combo:
                    for q, al in copies:
                        cnt[q][al] += 1
                    if trp >= 0:
                        tr[trp] += 1
                        kF[trp] += isibd
                    lc += c0
                groups[tuple(tuple(r) for r in cnt)][(tuple(kF), tuple(tr))] += math.exp(lc)
            fac = np.zeros((len(Fwt), len(Twt)))
            pim = np.zeros((Kl, len(Fwt), len(Twt)))
            pm = np.zeros((P, Kl, 2, len(Fwt), len(Twt)))
            for cnt, fk in groups.items():
                fsum = sum(w * fpoly(kF, tr) for (kF, tr), w in fk.items())
                like, pimom, pmom = dm_terms(li, cnt)
                fac += np.outer(fsum, like)
                pim += fsum[None, :, None] * pimom[:, None, :]
                pm += fsum[None, None, None, :, None] * pmom[:, :, :, None, :]
            loc.append((fac, pim, pm))

        prod = np.prod([x[0] for x in loc], axis=0)               # nF^P x nT^P
        W2 = np.outer(Fwt, Twt)
        like = float((W2 * prod).sum())
        if like <= 0:
            continue
        logw = logPA + math.log(like)
        if offset is None:
            offset = logw
        if logw - offset > 600:
            sc = math.exp(offset - logw)
            for k in acc:
                acc[k] *= sc
            offset = logw
        wt = math.exp(logw - offset)
        acc["Z"] += wt
        for (s, j), v in Em.items():
            acc[("m", s, j, 1)] += wt * v
            acc[("m", s, j, 2)] += wt * Em2[(s, j)]
            acc[("z0", s, j)] += wt * Z0[(s, j)]
        post = W2 * prod / like
        for q in range(P):
            fq = Fg[q][:, None]
            tq = Tg[q][None, :]
            acc[("F", q, 1)] += wt * float((post * fq).sum())
            acc[("F", q, 2)] += wt * float((post * fq ** 2).sum())
            acc[("theta", q, 1)] += wt * float((post * tq).sum())
            acc[("theta", q, 2)] += wt * float((post * tq ** 2).sum())
        for n in range(N):
            acc[("anc", n, anc[n])] += wt
        for li in range(L):
            fac, pim, pm = loc[li]
            others = W2 * prod / np.where(fac == 0, 1, fac) / like
            for k in range(K[li]):
                acc[("pi", li, k)] += wt * float((others * pim[k]).sum())
                for q in range(P):
                    for mom in (0, 1):
                        acc[("freq", q, li, k, mom + 1)] += wt * float((others * pm[q, k, mom]).sum())

    Z = acc["Z"]
    out = {}
    for s in range(P):
        for j in range(P):
            if s != j:
                m1, m2 = acc[("m", s, j, 1)] / Z, acc[("m", s, j, 2)] / Z
                out[f"m[{pname[s]}][{pname[j]}]"] = (m1, math.sqrt(max(m2 - m1 * m1, 0)))
                # Savage-Dickey BF_01 for m = 0: Pr(c = 0) (n_s+P-1)/(P-1)
                out[f"p0[{pname[s]}][{pname[j]}]"] = (acc[("z0", s, j)] / Z, None)
                out[f"bf[{pname[s]}][{pname[j]}]"] = (acc[("z0", s, j)] / Z * (samp.count(s) + P - 1) / (P - 1), None)
    for q in range(P):
        for nm in ("F", "theta"):
            f1, f2 = acc[(nm, q, 1)] / Z, acc[(nm, q, 2)] / Z
            out[f"{nm}[{pname[q]}]"] = (f1, math.sqrt(max(f2 - f1 * f1, 0)))
    for li, l in enumerate(loci):
        for k, name in enumerate(allele_list[l]):
            out[f"pi {l}/{name}"] = (acc[("pi", li, k)] / Z, None)
            for q in range(P):
                f1 = acc[("freq", q, li, k, 1)] / Z
                f2 = acc[("freq", q, li, k, 2)] / Z
                out[f"freq {pname[q]}/{l}/{name}"] = (f1, math.sqrt(max(f2 - f1 * f1, 0)))
    for n, i in enumerate(inds):
        for st in St[n]:
            out[f"anc {i}[{pname[st[0]]},{st[1]}]"] = (acc[("anc", n, st)] / Z, None)
    if a.json:
        json.dump({"post": out}, open(a.json, "w"), indent=1)
    for k in sorted(out):
        if k.startswith(("m[", "F[", "theta[", "bf[")):
            print(f"{k:14s} mean {out[k][0]:.5f}" + (f"  SD {out[k][1]:.5f}" if out[k][1] is not None else ""))
    print(f"({len(configs)} ancestry configurations; nodes F {nF}, theta {a.nodes_theta}, "
          f"pi {a.nodes_pi}/stick)")


if __name__ == "__main__":
    main()
