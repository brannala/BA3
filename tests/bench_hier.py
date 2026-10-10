#!/usr/bin/env python3
"""Speed and sampling efficiency of collapsed-sampler variants.

Runs each configuration (LABEL=BINARY[:FLAGS]) on each dataset with the same
seed and iterations, in parallel (one run per core; use an idle machine with
more cores than runs), and reports per configuration and dataset:

  sec/100k   wall-clock seconds per 100,000 iterations
  ESS        effective sample size of the off-diagonal migration-rate traces
             after burn-in (median and minimum over the m[i][j] columns),
             Geyer initial-positive-sequence autocorrelation time (bench.py)
  ESS/s      median ESS per second

usage: tests/bench_hier.py LABEL=BINARY[:FLAGS] ... --data NAME=ARGS ...
                           [--iters N] [--sampling K] [--seed S]
  e.g. tests/bench_hier.py mh=/tmp/350/BA3:-c gibbs=./BA3:-c "hier=./BA3:-c -H" \\
         --data "sc3=-V sc3.vcf -M sc3_meta.txt" --iters 500000
"""
import argparse, os, shutil, subprocess, sys, tempfile, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import iat  # noqa: E402


def run(binary, flags, args, iters, sampling, seed):
    d = tempfile.mkdtemp(prefix="ba3bench_")
    cmd = [binary] + flags + [os.path.abspath(a) if os.path.exists(a) else a for a in args] + \
          ["-s", str(seed), "-i", str(iters), "-b", str(iters // 10), "-n", str(sampling), "-t", "-o", "out.txt"]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    secs = time.time() - t0
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)}: {r.stderr[-400:]}")
    rows = [l.rstrip("\t\n").split("\t") for l in open(os.path.join(d, "BA3trace.txt"))]
    head, X = rows[0], np.array(rows[1:], float)
    X = X[X[:, 0] > iters // 10]
    cols = [k for k, h in enumerate(head) if h.startswith("m[") and h[2:].split("]")[0] != h.split("[")[2].rstrip("]")]
    ess = [len(X) / iat(X[:, k]) for k in cols]
    shutil.rmtree(d)
    return secs, float(np.median(ess)), float(np.min(ess))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("builds", nargs="+")
    ap.add_argument("--data", action="append", required=True)
    ap.add_argument("--iters", type=int, default=500_000)
    ap.add_argument("--sampling", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    builds = []
    for b in a.builds:
        label, rest = b.split("=", 1)
        binary, _, fl = rest.partition(":")
        builds.append((label, os.path.abspath(binary), fl.split()))
    data = [(d.split("=", 1)[0], d.split("=", 1)[1].split()) for d in a.data]
    jobs = [(dn, lab, bn, fl, args) for dn, args in data for lab, bn, fl in builds]
    with ThreadPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(lambda j: run(j[2], j[3], j[4], a.iters, a.sampling, a.seed), jobs))
    print(f"{'dataset':16} {'build':8} {'sec/100k':>9} {'ESS med':>8} {'ESS min':>8} {'ESS/s':>8}")
    for (dn, lab, *_), (secs, med, mn) in zip(jobs, res):
        print(f"{dn:16} {lab:8} {secs / a.iters * 1e5:9.2f} {med:8.0f} {mn:8.0f} {med / secs:8.2f}")


if __name__ == "__main__":
    main()
