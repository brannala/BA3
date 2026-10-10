#!/usr/bin/env python3
"""Prototype collapsed sampler for BA3 with a hierarchical allele-frequency prior.

Python reference for the C++ --hier work (doc/hierarchical_prior_plan.md).
Model as in tests/exact_posterior_hier.py; with --flat it is BA3's current model
(independent Dirichlet(1) frequencies per population), as in exact_posterior.py.

Allele frequencies are integrated out (Dirichlet-multinomial). State: each
individual's ancestry (source, age), the IBD indicator of each age-0/1
homozygote, the phase of each age-2 heterozygote, m, F, and (unless --flat) the
ancestral frequencies pi_l and drift theta_q, c_q = (1 - theta_q) / theta_q.

One sweep:
  * each individual in turn: remove its gene copies, draw its ancestry from the
    full conditional over all 2P-1 states (m-prior x collapsed predictive, IBD
    and phase marginalised), re-add it, drawing IBD indicators and phases;
  * m: Gibbs, each row of u = 3m is Dirichlet given the ancestry counts;
  * F: Gibbs, Beta given the IBD indicators;
  * pi_l: MH random walk in additive-log-ratio coordinates, all loci at once;
  * theta_q: MH random walk on logit(theta_q).

usage:
  tests/proto_hier.py data.txt [--flat] [--exact exact.json]      # BA3 format
  tests/proto_hier.py --vcf x.vcf --meta x_meta.txt [--flat]       # VCF (biallelic)
"""
import argparse, json, math, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_posterior import read_ba3  # noqa: E402

lgam = np.vectorize(math.lgamma, otypes=[float])


def load_native(fn):
    inds, pop_of, pops, loci, allele_sets, geno = read_ba3(fn)
    al = {l: sorted(allele_sets[l]) for l in loci}
    K = np.array([len(al[l]) for l in loci])
    G = -np.ones((len(inds), len(loci), 2), dtype=int)
    for n, i in enumerate(inds):
        for li, l in enumerate(loci):
            g = geno[(i, l)]
            if g is not None:
                G[n, li] = (al[l].index(g[0]), al[l].index(g[1]))
    samp = np.array([pops.index(pop_of[i]) for i in inds])
    return inds, pops, loci, al, K, G, samp


def load_vcf(vcf, meta):
    pop_of = dict(l.split()[:2] for l in open(meta) if l.strip())
    pops, rows, ids = [], [], None
    for line in open(vcf):
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            ids = f[9:]
            continue
        row = []
        for g in f[9:]:
            g = g.split(':')[0].replace('|', '/')
            row.append((-1, -1) if '.' in g else tuple(int(x) for x in g.split('/')))
        rows.append(row)
    keep = [k for k, i in enumerate(ids) if i in pop_of]
    inds = [ids[k] for k in keep]
    for i in inds:
        if pop_of[i] not in pops:
            pops.append(pop_of[i])
    G = np.array(rows, dtype=int)[:, keep, :].transpose(1, 0, 2)
    K = np.full(G.shape[1], 2)
    loci = [f"L{l}" for l in range(G.shape[1])]
    al = {l: ["0", "1"] for l in loci}
    samp = np.array([pops.index(pop_of[i]) for i in inds])
    return inds, pops, loci, al, K, G, samp


class Sampler:
    def __init__(self, G, samp, K, P, hier, rng):
        self.G, self.samp, self.K, self.P, self.hier, self.rng = G, samp, K, P, hier, rng
        N, L = G.shape[:2]
        self.N, self.L, self.A = N, L, int(K.max())
        self.valid = np.arange(self.A)[None, :] < K[:, None]            # L x A
        self.obs = G[:, :, 0] >= 0                                        # N x L
        self.hom = self.obs & (G[:, :, 0] == G[:, :, 1])
        self.het = self.obs & ~self.hom
        self.g0, self.g1 = np.maximum(G[:, :, 0], 0), np.maximum(G[:, :, 1], 0)
        self.Lidx = np.arange(L)
        self.src = samp.copy()
        self.age = np.zeros(N, dtype=int)
        self.z = np.zeros((N, L), dtype=bool)         # IBD (age 0/1 homozygotes)
        self.ph = np.zeros((N, L), dtype=bool)        # age-2 het: True = a1 from source
        self.cnt = np.zeros((P, L, self.A))
        self.m = np.full((P, P), 0.01)
        np.fill_diagonal(self.m, 0.0)
        self.F = np.full(P, 0.05)
        # pooled sample frequencies as the starting pi
        c = np.zeros((L, self.A))
        for k in (0, 1):
            np.add.at(c, (np.nonzero(self.obs)[1], G[:, :, k][self.obs]), 1)
        self.pi = np.where(self.valid, c + 0.5, 0.0)
        self.pi /= self.pi.sum(1, keepdims=True)
        self.theta = np.full(P, 0.1)
        self.sd_pi, self.sd_th = 0.5, 0.3
        self.acc = {'pi': [0, 0], 'theta': [0, 0]}
        for n in range(N):
            self._add(n)

    # ---- prior concentrations
    def alpha(self, q):
        """(L x A) Dirichlet parameters for population q, and their sums (L)."""
        if self.hier:
            c = (1 - self.theta[q]) / self.theta[q]
            return c * self.pi, np.full(self.L, c)
        a = self.valid.astype(float)
        return a, self.K.astype(float)

    # ---- count-table bookkeeping
    def _copies(self, n):
        """(pop, locus, allele) arrays of the gene copies individual n contributes."""
        P_, L_, A_ = [], [], []
        g0, g1, obs = self.g0[n], self.g1[n], self.obs[n]
        if self.age[n] < 2:
            q = self.src[n]
            one = obs & self.hom[n] & self.z[n]                  # IBD: one copy
            two = obs & ~(self.hom[n] & self.z[n])               # else two copies
            for msk, al in ((obs, g0), (two, g1)):
                L_.append(self.Lidx[msk]); A_.append(al[msk]); P_.append(np.full(msk.sum(), q))
        else:
            mig, nat = self.src[n], self.samp[n]
            sw = self.ph[n]
            a_src = np.where(sw, g1, g0)
            a_nat = np.where(sw, g0, g1)
            for q, al in ((mig, a_src), (nat, a_nat)):
                L_.append(self.Lidx[obs]); A_.append(al[obs]); P_.append(np.full(obs.sum(), q))
        return np.concatenate(P_), np.concatenate(L_), np.concatenate(A_)

    def _remove(self, n):
        p, l, a = self._copies(n)
        np.add.at(self.cnt, (p, l, a), -1)

    def _pred(self, n, src, age):
        """Log predictive of individual n's genotypes (n removed) under (src, age),
        IBD and phase marginalised, and the per-locus pieces needed to draw them."""
        g0, g1, obs, hom, het, L = self.g0[n], self.g1[n], self.obs[n], self.hom[n], self.het[n], self.Lidx
        if age < 2:
            q = src
            al, S = self.alpha(q)
            cn = self.cnt[q]
            Nq = cn.sum(1)
            n0, n1 = cn[L, g0] + al[L, g0], cn[L, g1] + al[L, g1]
            r = n0 / (Nq + S)
            F = self.F[q]
            pz1 = F * r
            pz0 = (1 - F) * r * (n0 + 1) / (Nq + 1 + S)
            vhom = pz1 + pz0
            vhet = (1 - F) * r * n1 / (Nq + 1 + S)
            v = np.where(hom, vhom, vhet)
            lp = np.log(v[obs]).sum()
            return lp, pz1 / np.where(hom, vhom, 1.0)
        mig, nat = src, self.samp[n]
        am, Sm = self.alpha(mig)
        an, Sn = self.alpha(nat)
        cm, cn_ = self.cnt[mig], self.cnt[nat]
        Dm, Dn = cm.sum(1) + Sm, cn_.sum(1) + Sn
        rm0, rm1 = (cm[L, g0] + am[L, g0]) / Dm, (cm[L, g1] + am[L, g1]) / Dm
        rn0, rn1 = (cn_[L, g0] + an[L, g0]) / Dn, (cn_[L, g1] + an[L, g1]) / Dn
        w0, w1 = rm0 * rn1, rm1 * rn0                    # w0: a0 from source
        v = np.where(hom, w0, 0.5 * (w0 + w1))
        lp = np.log(v[obs]).sum()
        return lp, w1 / np.where(het, w0 + w1, 1.0)

    def _add(self, n, extra=None):
        _, prob = extra if extra is not None else self._pred(n, self.src[n], self.age[n])
        u = self.rng.random(self.L)
        if self.age[n] < 2:
            self.z[n] = self.hom[n] & (u < prob)
            self.ph[n] = False
        else:
            self.ph[n] = self.het[n] & (u < prob)
            self.z[n] = False
        p, l, a = self._copies(n)
        np.add.at(self.cnt, (p, l, a), 1)

    # ---- moves
    def sweep_ancestry(self):
        P = self.P
        for n in self.rng.permutation(self.N):
            s = self.samp[n]
            self._remove(n)
            sts = [(s, 0)] + [(j, a) for j in range(P) if j != s for a in (1, 2)]
            pri = [1 - 3 * self.m[s].sum()] + [self.m[s, j] * (1 if a == 1 else 2) for j, a in sts[1:]]
            res = [self._pred(n, j, a) for j, a in sts]
            lw = np.array([math.log(pr) + r[0] if pr > 0 else -np.inf for pr, r in zip(pri, res)])
            w = np.exp(lw - lw.max())
            k = self.rng.choice(len(sts), p=w / w.sum())
            self.src[n], self.age[n] = sts[k]
            self._add(n, res[k])

    def gibbs_m(self):
        P = self.P
        for s in range(P):
            mem = self.samp == s
            others = [j for j in range(P) if j != s]
            cj = [np.sum(mem & (self.src == j) & (self.age > 0)) for j in others]
            n0 = np.sum(mem & (self.age == 0))
            u = self.rng.dirichlet(np.array(cj + [n0], dtype=float) + 1)
            for j, uj in zip(others, u[:-1]):
                self.m[s, j] = uj / 3

    def gibbs_F(self):
        for q in range(self.P):
            who = (self.age < 2) & (self.src == q)
            T = self.obs[who].sum()
            Z = self.z[who].sum()
            self.F[q] = self.rng.beta(1 + Z, 1 + T - Z)

    def _logdm(self, q, al, S):
        cn = self.cnt[q]
        t = np.where(self.valid, lgam(np.where(self.valid, al + cn, 1)) - lgam(np.where(self.valid, al, 1)), 0)
        return lgam(S) - lgam(S + cn.sum(1)) + t.sum(1)                 # per locus

    def mh_pi(self):
        if not self.hier:
            return
        y = np.where(self.valid, np.log(np.where(self.valid, self.pi, 1)), -np.inf)
        y = y - y[np.arange(self.L), self.K - 1][:, None]
        yp = y + np.where(self.valid, self.rng.normal(0, self.sd_pi, y.shape), 0)
        yp[np.arange(self.L), self.K - 1] = 0
        pp = np.exp(yp - np.nanmax(np.where(self.valid, yp, -np.inf), 1, keepdims=True))
        pp = np.where(self.valid, pp, 0)
        pp /= pp.sum(1, keepdims=True)

        def target(pi):
            t = np.zeros(self.L)
            for q in range(self.P):
                c = (1 - self.theta[q]) / self.theta[q]
                t += self._logdm(q, c * pi, np.full(self.L, c))
            return t + np.where(self.valid, np.log(np.where(self.valid, pi, 1)), 0).sum(1)   # alr Jacobian

        lr = target(pp) - target(self.pi)
        ok = np.log(self.rng.random(self.L)) < lr
        self.pi[ok] = pp[ok]
        self.acc['pi'][0] += ok.sum(); self.acc['pi'][1] += self.L

    def mh_theta(self):
        if not self.hier:
            return
        for q in range(self.P):
            th = self.theta[q]
            yp = math.log(th / (1 - th)) + self.rng.normal(0, self.sd_th)
            tp = 1 / (1 + math.exp(-yp))
            if not 1e-9 < tp < 1 - 1e-9:
                continue

            def target(t):
                c = (1 - t) / t
                return self._logdm(q, c * self.pi, np.full(self.L, c)).sum() + math.log(t * (1 - t))

            ok = math.log(self.rng.random()) < target(tp) - target(th)
            if ok:
                self.theta[q] = tp
            self.acc['theta'][0] += ok; self.acc['theta'][1] += 1

    def tune(self):
        for k, attr in (('pi', 'sd_pi'), ('theta', 'sd_th')):
            a, t = self.acc[k]
            if t:
                rate = a / t
                setattr(self, attr, getattr(self, attr) * (1.2 if rate > 0.4 else 0.8 if rate < 0.2 else 1.0))
            self.acc[k] = [0, 0]

    def sweep(self):
        self.sweep_ancestry()
        self.gibbs_m()
        self.gibbs_F()
        for _ in range(3):
            self.mh_pi()
        self.mh_theta()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data", nargs="?")
    ap.add_argument("--vcf"); ap.add_argument("--meta")
    ap.add_argument("--flat", action="store_true", help="BA3's current independent Dirichlet(1) prior")
    ap.add_argument("--sweeps", type=int, default=2000)
    ap.add_argument("--burn", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--exact", help="exact-posterior JSON to compare against")
    ap.add_argument("--json", default="")
    ap.add_argument("--every", type=int, default=0, help="print a progress line every k sweeps")
    a = ap.parse_args()

    if a.vcf:
        inds, pops, loci, al, K, G, samp = load_vcf(a.vcf, a.meta)
    else:
        inds, pops, loci, al, K, G, samp = load_native(a.data)
    P = len(pops)
    rng = np.random.default_rng(a.seed)
    S = Sampler(G, samp, K, P, not a.flat, rng)

    keys, rows = None, []
    for it in range(1, a.sweeps + 1):
        S.sweep()
        if it <= a.burn and it % 50 == 0:
            S.tune()
        if a.every and it % a.every == 0:
            print(f"sweep {it}: age2 {np.bincount(S.age, minlength=3).tolist()} "
                  f"maxm {S.m.max():.3f} F {np.round(S.F, 3).tolist()} theta {np.round(S.theta, 3).tolist()}",
                  flush=True)
        if it <= a.burn:
            continue
        row = {}
        for s in range(P):
            for j in range(P):
                if s != j:
                    row[f"m[{pops[s]}][{pops[j]}]"] = S.m[s, j]
        for q in range(P):
            row[f"F[{pops[q]}]"] = S.F[q]
            if S.hier:
                row[f"theta[{pops[q]}]"] = S.theta[q]
        if S.N * S.L <= 2000:          # per-locus and per-individual output on small data only
            for q in range(P):
                alq, Sq = S.alpha(q)
                ef = (S.cnt[q] + alq) / (S.cnt[q].sum(1) + Sq)[:, None]
                for li, l in enumerate(loci):
                    for k, nm in enumerate(al[l]):
                        row[f"freq {pops[q]}/{l}/{nm}"] = ef[li, k]
            for n, i in enumerate(inds):
                s = samp[n]
                for j, ag in [(s, 0)] + [(j, x) for j in range(P) if j != s for x in (1, 2)]:
                    row[f"anc {i}[{pops[j]},{ag}]"] = float(S.src[n] == j and S.age[n] == ag)
        keys = keys or list(row)
        rows.append([row[k] for k in keys])

    X = np.array(rows)
    mean, sd = X.mean(0), X.std(0)
    nb = 20
    bm = np.array([b.mean(0) for b in np.array_split(X, nb)])
    se = bm.std(0, ddof=1) / math.sqrt(nb)                     # batch-means SE
    post = {k: (float(mu), float(s), float(e)) for k, mu, s, e in zip(keys, mean, sd, se)}
    if a.json:
        json.dump({"post": post}, open(a.json, "w"), indent=1)
    for k in keys:
        if k.startswith(("m[", "F[", "theta[")):
            print(f"{k:16s} mean {post[k][0]:.4f}  SD {post[k][1]:.4f}  (SE {post[k][2]:.4f})")
    print(f"ages 0/1/2 at end: {np.bincount(S.age, minlength=3).tolist()}")
    if a.exact:
        ex = json.load(open(a.exact))["post"]
        worst = []
        for k in keys:
            if k in ex and ex[k][0] is not None:
                z = (post[k][0] - ex[k][0]) / max(post[k][2], 1e-4)
                worst.append((abs(z), k, post[k][0], ex[k][0]))
        worst.sort(reverse=True)
        print(f"vs exact: {len(worst)} quantities; max |z| {worst[0][0]:.1f}; "
              f"{sum(w[0] > 3 for w in worst)} with |z| > 3")
        for w in worst[:6]:
            print(f"   |z| {w[0]:4.1f}  {w[1]}: proto {w[2]:.4f}  exact {w[3]:.4f}")


if __name__ == "__main__":
    main()
