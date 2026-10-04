#!/usr/bin/env python3
"""Speed and sampling-efficiency benchmark: standard sampler vs --collapse.

For each dataset, runs the standard sampler (the reference binary, normally the
master build) and the collapsed sampler (the new binary with -c) once each, with
the same seed and iteration count, one run at a time so the timings don't
compete for cores. Reports:

  sec/100k   wall-clock seconds per 100,000 MCMC iterations
  speedup    standard time / collapse time (per iteration)
  ESS        effective sample size of the off-diagonal migration-rate traces
             after burn-in (median and minimum over the m[i][j] columns),
             via the Geyer initial-positive-sequence autocorrelation time
  ESS/s      median ESS per second of wall-clock time -- the end-to-end
             efficiency measure (it combines per-iteration cost and mixing)

requires numpy.  usage:
  tests/bench.py <standard BA3> <collapse BA3> [--only NAME,...] [--scale X]
                 [--snp DIR]      (DIR = examples/snp, adds the empirical SNP sets)
"""
import argparse, os, shutil, subprocess, sys, tempfile, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX, TD = os.path.join(ROOT, "examples"), os.path.join(ROOT, "tests", "data")

# name, input args, iterations (burn-in = 10%)
DATASETS = [
    ("2pop (300 ind, 10 loci)",          [f"{EX}/2pop.txt"],                        2_000_000),
    ("3pop (400 ind, 7 loci)",           [f"{EX}/3pop.txt"],                        2_000_000),
    ("100loci_mig (200 ind, 100 loci)",  [f"{EX}/100loci_with_migration.txt"],      1_000_000),
    ("100loci_missing",                  [f"{TD}/100loci_missing.txt"],             1_000_000),
    ("very_high_info (60 ind, 20 loci)", [f"{EX}/very_high_info_no_migration.txt"], 2_000_000),
    ("litt VCF (106 ind, 12718 SNPs)",   ["-V", f"{EX}/allpopslitt.vcf", "-M", f"{EX}/allpopslitt_meta.txt"],
                                                                                      200_000),
]
SNP = [  # relative to --snp DIR
    ("dace (145 ind, 14355 SNPs)", "dace/spd_dv.filt2.vcf", "dace/spd_dv_meta.txt", 200_000),
    ("palm Bm (27 ind, 1935 SNPs)", "palm/Bm_miss50.vcf", "palm/Bm_meta.txt", 500_000),
    ("palm HC (43 ind, 12913 SNPs)", "palm/HC_miss50.vcf", "palm/HC_meta.txt", 200_000),
]


def iat(x):
    """Integrated autocorrelation time, Geyer initial positive sequence."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    if n < 4 or np.dot(x, x) == 0:
        return float("nan")
    f = np.fft.rfft(x, n=2 * n)
    acov = np.fft.irfft(f * np.conj(f))[:n].real
    rho = acov / acov[0]
    tau, t = 1.0, 1
    while t + 1 < n:
        pair = rho[t] + rho[t + 1]
        if pair <= 0:
            break
        tau += 2.0 * pair
        t += 2
    return max(tau, 1.0)


def run(binary, collapse, args, iters, thin, seed):
    d = tempfile.mkdtemp(prefix="ba3bench_")
    cmd = [binary] + (["-c"] if collapse else []) + [
        "-s", str(seed), "-i", str(iters), "-b", str(iters // 10), "-n", str(thin),
        "-t", "-o", "out.txt"] + args
    t0 = time.perf_counter()
    res = subprocess.run(cmd, cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    wall = time.perf_counter() - t0
    if res.returncode != 0:
        sys.exit(f"{' '.join(cmd)} failed:\n{res.stderr[-2000:]}")
    fn = os.path.join(d, "BA3trace.txt")
    hdr = open(fn).readline().rstrip("\n").split("\t")
    cols = [i for i, h in enumerate(hdr)
            if h.startswith("m[") and h.split("][")[0][2:] != h.split("][")[1][:-1]]
    data = np.loadtxt(fn, skiprows=1, usecols=cols, ndmin=2)
    shutil.rmtree(d)
    data = data[len(data) // 10:]              # drop burn-in samples
    ess = [len(c) / iat(c) for c in data.T if np.std(c) > 0]
    return wall, np.array(ess)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("std_bin"); ap.add_argument("col_bin")
    ap.add_argument("--only", default="")
    ap.add_argument("--scale", type=float, default=1.0, help="multiply iteration counts")
    ap.add_argument("--snp", default="", help="examples/snp directory (adds the SNP sets)")
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    sets = list(DATASETS)
    if a.snp:
        sets += [(n, ["-V", os.path.join(a.snp, v), "-M", os.path.join(a.snp, m)], it)
                 for n, v, m, it in SNP]
    only = [s for s in a.only.split(",") if s]
    std_bin, col_bin = os.path.abspath(a.std_bin), os.path.abspath(a.col_bin)

    print(f"{'dataset':<34} {'mode':<9} {'iters':>9} {'sec/100k':>9} {'speedup':>8} "
          f"{'ESS med':>8} {'ESS min':>8} {'ESS/s':>8} {'ESS/s x':>8}")
    for name, args, iters in sets:
        if only and not any(o in name for o in only):
            continue
        iters = int(iters * a.scale)
        thin = max(1, iters // 20000)          # ~20k trace samples per run
        rows = {}
        for mode, b, c in (("standard", std_bin, False), ("collapse", col_bin, True)):
            wall, ess = run(b, c, args, iters, thin, a.seed)
            rows[mode] = (wall, ess)
        ws, es = rows["standard"]; wc, ec = rows["collapse"]
        for mode in ("standard", "collapse"):
            w, e = rows[mode]
            med = np.median(e) if len(e) else float("nan")
            mn = e.min() if len(e) else float("nan")
            sp = f"{ws / wc:8.1f}" if mode == "collapse" else f"{'':>8}"
            eff = (med / w) / (np.median(es) / ws) if mode == "collapse" and len(es) else None
            ex = f"{eff:8.1f}" if eff is not None else f"{'':>8}"
            print(f"{name if mode == 'standard' else '':<34} {mode:<9} {iters:>9} "
                  f"{w / iters * 1e5:9.2f} {sp} {med:8.0f} {mn:8.0f} {med / w:8.1f} {ex}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
