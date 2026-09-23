#!/usr/bin/env python3
"""Posterior-agreement test: --collapse vs the standard sampler.

The collapsed sampler targets the same posterior as the standard sampler (allele
frequencies are integrated out analytically instead of sampled), but it is a
different Markov chain, so the two cannot be compared seed-for-seed. Instead each
mode is run with R independent seeds and every reported posterior mean is compared
using the between-seed Monte Carlo standard error:

    z = (mean_collapse - mean_standard) / sqrt(se_c^2 + se_s^2 + floor^2)

where se is the SD of the per-seed posterior means / sqrt(R) and `floor` absorbs
the output rounding (e.g. ancestry probabilities are printed to 3 decimals).

Quantities compared: migration rates m[i][j], inbreeding F, allele frequencies
(-F file), and per-individual migrant-ancestry probabilities (BA3indiv.txt).

If both samplers target the same posterior, z is approximately t-distributed, so
~5% of |z| > 2 and very few > 3. A group FAILS if more than 10% of its |z| exceed
2 or more than 2% exceed 3.

usage: tests/compare_collapse.py <standard BA3> <collapse BA3> [--seeds R]
         [--jobs J] [--quick] [--only NAME[,NAME...]]
The <standard BA3> is normally the master build (the reference implementation);
the <collapse BA3> is run with -c.
"""
import argparse, math, os, re, shutil, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX, TD = os.path.join(ROOT, "examples"), os.path.join(ROOT, "tests", "data")

# name, input args, iterations, burn-in (quick mode divides iterations by 10)
DATASETS = [
    ("2pop",            [f"{EX}/2pop.txt"],                         1_000_000, 100_000),
    ("3pop",            [f"{EX}/3pop.txt"],                         1_000_000, 100_000),
    ("100loci_mig",     [f"{EX}/100loci_with_migration.txt"],       1_000_000, 100_000),
    ("high_info",       [f"{EX}/high_info_no_migration.txt"],       1_000_000, 100_000),
    ("very_high_info",  [f"{EX}/very_high_info_no_migration.txt"],  1_000_000, 100_000),
    ("2pop_missing",    [f"{TD}/2pop_missing.txt"],                 1_000_000, 100_000),
    ("100loci_missing", [f"{TD}/100loci_missing.txt"],              1_000_000, 100_000),
    ("vcf_litt",        ["-V", f"{EX}/allpopslitt.vcf", "-M", f"{EX}/allpopslitt_meta.txt"],
                                                                      500_000,  50_000),
]
FLOOR = {"m": 5e-4, "F": 5e-4, "freq": 1e-5, "anc": 5e-4}


def parse(run_dir):
    """Return {(group, key): posterior mean} for one run."""
    out = {}
    txt = open(os.path.join(run_dir, "out.txt")).read().splitlines()
    for k, line in enumerate(txt):
        if "Migration Rate Matrix" in line:
            r = k + 3
            while r < len(txt) and re.match(r"\s*\[\d+\]", txt[r]):
                i = int(re.search(r"\[(\d+)\]", txt[r]).group(1))
                for j, (m, _) in enumerate(re.findall(r"([0-9.]+)\(([0-9.]+)\)", txt[r])):
                    if i != j:
                        out[("m", f"m[{i}][{j}]")] = float(m)
                r += 1
        if "Inbreeding Coefficients" in line:
            r = k + 3
            while r < len(txt) and re.match(r"\s*\[\s*\d+\]", txt[r]):
                i = int(re.search(r"\[\s*(\d+)\]", txt[r]).group(1))
                out[("F", f"F[{i}]")] = float(re.search(r"([0-9.]+)\(", txt[r]).group(1))
                r += 1
    with open(os.path.join(run_dir, "freq.txt")) as f:
        next(f)
        for line in f:
            p, l, a, fr, _ = line.split("\t")
            out[("freq", f"{p}/{l}/{a}")] = float(fr)
    ind = None
    for line in open(os.path.join(run_dir, "BA3indiv.txt")):
        m = re.match(r"\s*Individual: (\S+)", line)
        if m:
            ind = m.group(1)
            continue
        for pop, age, pr in re.findall(r"\[(\d+),(\d+)\]:([0-9.]+)", line):
            out[("anc", f"{ind}[{pop},{age}]")] = float(pr)
    return out


def run_one(binary, collapse, args, iters, burn, seed, work):
    d = tempfile.mkdtemp(dir=work)
    cmd = [binary] + (["-c"] if collapse else []) + [
        "-s", str(seed), "-i", str(iters), "-b", str(burn), "-n", "100",
        "-g", "-F", "freq.txt", "-o", "out.txt"] + args
    res = subprocess.run(cmd, cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed:\n{res.stderr[-2000:]}")
    vals = parse(d)
    shutil.rmtree(d)
    return vals


def summarize(runs):
    keys = runs[0].keys()
    R = len(runs)
    stats = {}
    for k in keys:
        xs = [r[k] for r in runs]
        mu = sum(xs) / R
        sd = math.sqrt(sum((x - mu) ** 2 for x in xs) / (R - 1)) if R > 1 else 0.0
        stats[k] = (mu, sd / math.sqrt(R))
    return stats


def compare(std, col):
    groups = {}
    for k in std:
        if k not in col:
            continue
        (ms, ss), (mc, sc) = std[k], col[k]
        fl = FLOOR[k[0]]
        z = (mc - ms) / math.sqrt(ss * ss + sc * sc + fl * fl)
        groups.setdefault(k[0], []).append((abs(z), abs(mc - ms), k[1], ms, mc))
    return groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("std_bin"); ap.add_argument("col_bin")
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--quick", action="store_true", help="iterations / 10 (smoke test)")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    std_bin, col_bin = os.path.abspath(a.std_bin), os.path.abspath(a.col_bin)
    only = set(a.only.split(",")) - {""}
    work = tempfile.mkdtemp(prefix="ba3cmp_")
    allok = True
    try:
        for name, args, iters, burn in DATASETS:
            if only and name not in only:
                continue
            if a.quick:
                iters, burn = iters // 10, burn // 10
            with ThreadPoolExecutor(a.jobs) as ex:
                fs = [ex.submit(run_one, b, c, args, iters, burn, 1000 + s, work)
                      for (b, c) in ((std_bin, False), (col_bin, True)) for s in range(a.seeds)]
                res = [f.result() for f in fs]
            std, col = summarize(res[:a.seeds]), summarize(res[a.seeds:])
            print(f"\n== {name}  ({a.seeds} seeds x {iters} iters per mode)")
            for g, zs in sorted(compare(std, col).items()):
                zs.sort(reverse=True)
                n = len(zs)
                f2 = sum(z > 2 for z, *_ in zs) / n
                f3 = sum(z > 3 for z, *_ in zs) / n
                ok = f2 <= 0.10 and f3 <= 0.02
                allok &= ok
                zmax, dmax, key, ms, mc = zs[0]
                print(f"  {'PASS' if ok else 'FAIL'}  {g:<5} n={n:<7} |z|>2: {100*f2:5.1f}%  "
                      f"|z|>3: {100*f3:5.2f}%  max|z|={zmax:5.2f} ({key}: std {ms:.4f} col {mc:.4f})  "
                      f"max|diff|={max(d for _, d, *_ in zs):.4f}")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("\ncompare_collapse:", "PASS" if allok else "FAIL")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
