#!/usr/bin/env python3
"""Check BA3 builds against the exact posterior of tiny data sets.

For each data set in tests/data/exact_*.txt, computes the exact posterior with
tests/exact_posterior.py, runs every build R times (different seeds), and compares
each reported posterior mean and SD -- migration rates, F, allele frequencies
(-F file) and individual ancestry probabilities (BA3indiv.txt) -- with the exact
value:  z = (mean over seeds - exact) / sqrt(SE^2 + floor^2), where SE is the
between-seed standard error and floor absorbs the printed rounding. A group fails
if any |z| > 4 or more than 2% of its |z| exceed 3.

usage: tests/check_exact.py LABEL=BINARY[:FLAGS] ... [--seeds R] [--iters N]
                            [--jobs J] [--only NAME,...]
  e.g. tests/check_exact.py master=/tmp/ref/BA3 new=./BA3 collapse=./BA3:-c
"""
import argparse, glob, json, math, os, re, shutil, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
FLOOR = {"m": 5e-5, "F": 5e-5, "freq": 5e-6, "anc": 5e-4}


def parse(d):
    """Posterior (mean, sd) by name-based key, from one BA3 run directory."""
    txt = open(os.path.join(d, "out.txt")).read()
    labels = dict(re.findall(r"^\s*\[(\d+)\]\s+(\S+)\s*$", txt.split("Migration Rate Matrix")[0], re.M))
    out = {}
    mat = txt.split("Migration Rate Matrix")[1].split("Savage")[0].splitlines()
    for line in mat:
        mm = re.match(r"\s*\[(\d+)\]\s+(.*)", line)
        if not mm or "(" not in line:
            continue
        i = mm.group(1)
        for j, (m, s) in enumerate(re.findall(r"([0-9.]+)\(([0-9.]+)\)", mm.group(2))):
            if str(j) != i:
                out[f"m[{labels[i]}][{labels[str(j)]}]"] = (float(m), float(s))
    for q, name, m, s in re.findall(r"\[\s*(\d+)\]\s+(\S+)\s+([0-9.]+)\(([0-9.]+)\)", txt.split("Inbreeding Coefficients")[1]):
        out[f"F[{name}]"] = (float(m), float(s))
    with open(os.path.join(d, "freq.txt")) as f:
        next(f)
        for line in f:
            p, l, a, fr, sd = line.rstrip("\n").split("\t")
            out[f"freq {p}/{l}/{a}"] = (float(fr), float(sd))
    ind = None
    for line in open(os.path.join(d, "BA3indiv.txt")):
        mm = re.match(r"\s*Individual: (\S+)", line)
        if mm:
            ind = mm.group(1)
            continue
        for pop, age, pr in re.findall(r"\[(\d+),(\d+)\]:([0-9.]+)", line):
            out[f"anc {ind}[{labels[pop]},{age}]"] = (float(pr), None)
    return out


def run(binary, flags, data, seed, iters, work):
    d = tempfile.mkdtemp(dir=work)
    shutil.copy(data, os.path.join(d, "in.txt"))
    cmd = [binary] + flags + ["-s", str(seed), "-i", str(iters), "-b", str(iters // 10),
                              "-n", "20", "-g", "-F", "freq.txt", "-o", "out.txt", "in.txt"]
    r = subprocess.run(cmd, cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)}: {r.stderr[-500:]}")
    v = parse(d)
    shutil.rmtree(d)
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("builds", nargs="+")
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--iters", type=int, default=20_000_000)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--only", default="")
    ap.add_argument("--python", default=sys.executable, help="python with numpy for the enumerator")
    a = ap.parse_args()
    builds = []
    for b in a.builds:
        label, rest = b.split("=", 1)
        binary, _, fl = rest.partition(":")
        builds.append((label, os.path.abspath(binary), fl.split() if fl else []))
    only = set(a.only.split(",")) - {""}
    work = tempfile.mkdtemp(prefix="ba3exact_")
    allok = True
    try:
        for data in sorted(glob.glob(os.path.join(HERE, "data", "exact_*.txt"))):
            name = os.path.basename(data)[:-4]
            if only and name not in only:
                continue
            ej = os.path.join(work, name + ".json")
            subprocess.run([a.python, os.path.join(HERE, "exact_posterior.py"), data, "--json", ej],
                           check=True, stdout=subprocess.DEVNULL)
            exact = {k: tuple(v) for k, v in json.load(open(ej))["post"].items()}
            with ThreadPoolExecutor(a.jobs) as ex:
                futs = {(lab, s): ex.submit(run, bn, fl, data, 1 + s, a.iters, work)
                        for lab, bn, fl in builds for s in range(a.seeds)}
                res = {k: f.result() for k, f in futs.items()}
            print(f"\n== {name}  ({a.seeds} seeds x {a.iters:,} iterations)")
            for lab, _, _ in builds:
                runs = [res[(lab, s)] for s in range(a.seeds)]
                groups = {}
                for key, (emean, esd) in exact.items():
                    if key not in runs[0]:
                        continue
                    g = key.split()[0].split("[")[0]
                    for stat, ev in (("mean", emean), ("sd", esd)):
                        if ev is None:
                            continue
                        xs = [r[key][0 if stat == "mean" else 1] for r in runs]
                        mu = sum(xs) / len(xs)
                        se = math.sqrt(sum((x - mu) ** 2 for x in xs) / (len(xs) - 1) / len(xs))
                        z = (mu - ev) / math.sqrt(se * se + FLOOR[g] ** 2)
                        groups.setdefault((g, stat), []).append((abs(z), z, key, mu, ev))
                line = []
                for (g, stat), zs in sorted(groups.items()):
                    zs.sort(reverse=True)
                    f3 = sum(z > 3 for z, *_ in zs) / len(zs)
                    ok = zs[0][0] <= 4 and f3 <= 0.02
                    allok &= ok
                    line.append(f"{g}.{stat} {'ok ' if ok else 'BAD'} max|z|={zs[0][0]:5.1f}")
                    if not ok:
                        _, z, key, mu, ev = zs[0]
                        line[-1] += f" ({key}: {mu:.4f} vs exact {ev:.4f})"
                print(f"  {lab:10s} " + " | ".join(line))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("\ncheck_exact:", "PASS" if allok else "FAIL")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
