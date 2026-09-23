#!/usr/bin/env python3
"""Check BA3's model-comparison estimators against exact values.

  --ss      stepping-stone log marginal likelihood, for each tests/data/exact_*.txt
            and for exact_3pop with popB and popC pooled; also the pooled-vs-full
            log Bayes factor.
  --sdpool  Savage-Dickey BF01 for H0: equal allele frequencies of a population pair.

Exact values come from tests/exact_posterior.py (log evidence on BA3's ordered-copy
scale; --pool and --tie for the pooled and nested models). Each estimate is the
mean over R seeds, z = (mean - exact) / SE.

usage: tests/check_evidence.py BA3 [--seeds R] [--ss-iters N] [--sd-iters N] [--jobs J]
"""
import argparse, glob, json, math, os, re, shutil, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
# An evidence estimate also passes if within TOL nats of exact: stepping-stone with
# 16 rungs has a small systematic bias (a few 1e-3 nats here), irrelevant for
# model choice, that long runs can resolve statistically.
TOL = 0.01
SD_CASES = [("exact_3pop", "popB", "popC"), ("exact_3pop", "popA", "popB"), ("exact_k234", "pop0", "pop1")]


def exact(py, data, extra):
    with tempfile.NamedTemporaryFile(suffix=".json") as f:
        subprocess.run([py, os.path.join(HERE, "exact_posterior.py"), data, "--json", f.name] + extra,
                       check=True, stdout=subprocess.DEVNULL)
        return json.load(open(f.name))["logZ_ba3"]


def pooled_copy(data, a, b, work):
    out = os.path.join(work, os.path.basename(data)[:-4] + f"_pool_{a}_{b}.txt")
    with open(data) as f, open(out, "w") as g:
        for line in f:
            t = line.split()
            if len(t) == 5 and t[1] == b:
                t[1] = a
            g.write("   ".join(t) + "\n")
    return out


def run(ba3, flag, data, seed, iters, work):
    d = tempfile.mkdtemp(dir=work)
    shutil.copy(data, os.path.join(d, "in.txt"))
    r = subprocess.run([ba3, flag, "-s", str(seed), "-i", str(iters), "-b", str(iters // 10), "-n", "10",
                        "-o", "out.txt", "in.txt"], cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr[-500:])
    txt = open(os.path.join(d, "out.txt")).read()
    shutil.rmtree(d)
    return txt


def mse(xs):
    m = sum(xs) / len(xs)
    return m, math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1) / len(xs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ba3")
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--ss-iters", type=int, default=8_000_000)
    ap.add_argument("--sd-iters", type=int, default=4_000_000)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()
    ba3 = os.path.abspath(a.ba3)
    work = tempfile.mkdtemp(prefix="ba3ev_")
    ok = True
    try:
        # ---- stepping-stone
        cases = [(os.path.basename(f)[:-4], f, []) for f in sorted(glob.glob(os.path.join(DATA, "exact_*.txt")))]
        cases.append(("exact_3pop pooled B+C", pooled_copy(os.path.join(DATA, "exact_3pop.txt"), "popB", "popC", work), []))
        with ThreadPoolExecutor(a.jobs) as ex:
            ex_z = {name: ex.submit(exact, a.python, f, extra) for name, f, extra in cases}
            runs = {(name, s): ex.submit(run, ba3, "-S", f, 1 + s, a.ss_iters, work)
                    for name, f, _ in cases for s in range(a.seeds)}
            ex_z = {k: v.result() for k, v in ex_z.items()}
            est = {}
            for name, _, _ in cases:
                xs = [float(re.search(r"Stepping-stone log marginal likelihood = (-?[0-9.]+)", runs[(name, s)].result()).group(1))
                      for s in range(a.seeds)]
                est[name] = xs
        print(f"== stepping-stone log p(G)  ({a.seeds} seeds x {a.ss_iters:,} iterations)")
        for name, _, _ in cases:
            m, se = mse(est[name]); z = (m - ex_z[name]) / se if se > 0 else float("inf")
            good = abs(z) < 4 or abs(m - ex_z[name]) < TOL
            ok &= good
            print(f"  {'ok ' if good else 'BAD'} {name:24s} exact {ex_z[name]:10.4f}  BA3 {m:10.4f} +- {se:.4f}  diff {m - ex_z[name]:+.4f}  z {z:+.1f}")
        full, pool = "exact_3pop", "exact_3pop pooled B+C"
        d = [p - f for p, f in zip(est[pool], est[full])]
        m, se = mse(d); ex_bf = ex_z[pool] - ex_z[full]
        print(f"  log BF pooled B+C vs full: exact {ex_bf:+.4f}  BA3 {m:+.4f} +- {se:.4f}  z {(m - ex_bf) / se:+.1f}")
        ok &= abs((m - ex_bf) / se) < 4

        # ---- Savage-Dickey pooling test
        with ThreadPoolExecutor(a.jobs) as ex:
            exs = {c: ex.submit(lambda c=c: (exact(a.python, os.path.join(DATA, c[0] + ".txt"), ["--tie", f"{c[1]},{c[2]}"])
                                            - exact(a.python, os.path.join(DATA, c[0] + ".txt"), [])) / math.log(10))
                   for c in SD_CASES}
            sruns = {(ds, s): ex.submit(run, ba3, "-P", os.path.join(DATA, ds + ".txt"), 1 + s, a.sd_iters, work)
                     for ds in {c[0] for c in SD_CASES} for s in range(a.seeds)}
            exs = {k: v.result() for k, v in exs.items()}
        print(f"\n== Savage-Dickey log10 BF01, H0: equal allele frequencies  ({a.seeds} seeds x {a.sd_iters:,} iterations)")
        for c in SD_CASES:
            ds, pa, pb = c
            xs, sds = [], []
            for s in range(a.seeds):
                txt = sruns[(ds, s)].result()
                mm = re.search(r"\] (\S+) - \[\d+\] (\S+)\s+(-?[0-9.]+)\s+([0-9.]+)", "")
                for line in txt.split("Pooling test")[1].splitlines():
                    mm = re.match(r"\s*\[\d+\] (\S+) - \[\d+\] (\S+)\s+(-?[0-9.]+)\s+([0-9.]+)", line)
                    if mm and {mm.group(1), mm.group(2)} == {pa, pb}:
                        xs.append(float(mm.group(3))); sds.append(float(mm.group(4)))
            m, se = mse(xs); z = (m - exs[c]) / math.sqrt(se * se + 0.0005 ** 2)
            good = abs(z) < 4
            ok &= good
            print(f"  {'ok ' if good else 'BAD'} {ds} {pa}-{pb}: exact {exs[c]:+.4f}  BA3 {m:+.4f} +- {se:.4f}  z {z:+.1f}"
                  f"   (SD of ln ratio {sum(sds) / len(sds):.2f})")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("\ncheck_evidence:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
