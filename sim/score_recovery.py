#!/usr/bin/env python3
"""
score_recovery.py

Score BA3's estimates against the ground truth produced by simulate_sexbias.py.

Reads:
  <prefix>_truth.json   ground truth (from simulate_sexbias.py)
  BA3out.txt            BA3's result file (default; override with --ba3out)

Compares:
  * Migration rates: BA3's off-diagonal m[i][j] (fraction of pop i that are
    first-generation migrants from j) against the realized age-1 migrant
    fraction in the truth file.
  * Dispersal sex bias: BA3's estimated phi (female fraction of migrants), if
    present in the output, against the true phi and the empirical phi. BA3 must
    print a line of the form
        Dispersal female fraction phi: 0.7912(0.0345)
    (any line containing "phi" with a NUMBER(SD) token is accepted).

Prints a report and exits 0 if every compared quantity lies within
--tol standard deviations of its target, else 1 (handy for automated tests).

Usage:
  python3 sim/score_recovery.py --truth sim_run_truth.json --ba3out BA3out.txt
"""

import argparse
import json
import re
import sys

# NUMBER(SD) token, e.g. 0.0250(0.0123)
EST_RE = re.compile(r"([-+]?\d*\.?\d+)\s*\(\s*([-+]?\d*\.?\d+)\s*\)")
ROW_RE = re.compile(r"^\s*\[(\d+)\]")


def parse_migration_matrix(lines):
    """Parse BA3's 'Migration Rate Matrix m[i][j]' block into
    matrix[i][j] = (mean, sd). Returns {} if not found."""
    start = None
    for idx, line in enumerate(lines):
        if "Migration Rate Matrix" in line:
            start = idx
            break
    if start is None:
        return {}
    matrix = {}
    for line in lines[start + 1:]:
        m = ROW_RE.match(line)
        if not m:
            # stop once we leave the contiguous block of [i] rows
            if matrix:
                break
            continue
        i = int(m.group(1))
        ests = EST_RE.findall(line)
        matrix[i] = [(float(mean), float(sd)) for (mean, sd) in ests]
    return matrix


def parse_sexbias(lines):
    """Parse the sex-biased dispersal block. Each parameter line begins with the
    parameter name (rho / phiMove / phiBreed). Returns {name: (mean, sd)}."""
    out = {}
    for line in lines:
        toks = line.strip().split()
        if not toks:
            continue
        key = toks[0]
        if key in ("rho", "phiMove", "phiBreed") and key not in out:
            m = EST_RE.search(line)
            if m:
                out[key] = (float(m.group(1)), float(m.group(2)))
    return out


def realized_migration(truth):
    """Realized first-generation migrant fractions from the truth file:
    m_real[i][j] = (# age-1 individuals in pop i sourced from j) / (n in pop i).
    Also returns per-pop native (age-0) fractions."""
    npop = truth["params"]["npop"]
    n_by_pop = [0] * npop
    age1 = [[0] * npop for _ in range(npop)]
    native = [0] * npop
    for ind in truth["individuals"]:
        i = ind["pop_index"]
        n_by_pop[i] += 1
        if ind["age"] == 0:
            native[i] += 1
        elif ind["age"] == 1:
            age1[i][ind["source_index"]] += 1
    m_real = [[(age1[i][j] / n_by_pop[i] if n_by_pop[i] else 0.0)
               for j in range(npop)] for i in range(npop)]
    native_frac = [native[i] / n_by_pop[i] if n_by_pop[i] else 0.0
                   for i in range(npop)]
    return m_real, native_frac


def main():
    p = argparse.ArgumentParser(description="Score BA3 output vs simulate_sexbias.py ground truth.")
    p.add_argument("--truth", required=True, help="<prefix>_truth.json from the simulator")
    p.add_argument("--ba3out", default="BA3out.txt", help="BA3 result file")
    p.add_argument("--tol", type=float, default=2.0,
                   help="pass if target within this many SDs of the estimate")
    p.add_argument("--json", action="store_true", help="also emit a JSON report to stdout")
    args = p.parse_args()

    truth = json.load(open(args.truth))
    try:
        lines = open(args.ba3out).read().splitlines()
    except OSError as e:
        sys.exit("Error opening BA3 output '%s': %s" % (args.ba3out, e))

    matrix = parse_migration_matrix(lines)
    sexbias = parse_sexbias(lines)
    m_real, native_frac = realized_migration(truth)
    npop = truth["params"]["npop"]

    report = {"migration": [], "sexbias": {}, "pass": True}

    print("=" * 64)
    print("BA3 recovery vs ground truth")
    print("=" * 64)

    # ---- migration rates (off-diagonal) ----
    true_m = truth["params"].get("m")
    if true_m is not None:
        print("\nMigration rate m[i][j], i<-j  (true m = %.4f per source; "
              "age-2 rate is 2m):" % true_m)
    else:
        print("\nMigration rates  m[i][j], i<-j  (off-diagonal only):")
    print("  %-8s %-14s %-14s %-8s %s" % ("i<-j", "BA3 mean(SD)", "realized age-1", "z", "within tol"))
    if not matrix:
        print("  (no migration matrix found in %s)" % args.ba3out)
    for i in range(npop):
        for j in range(npop):
            if i == j:
                continue
            realized = m_real[i][j]
            est = matrix.get(i, [])
            if j < len(est):
                mean, sd = est[j]
                z = (realized - mean) / sd if sd > 0 else float("inf")
                ok = abs(z) <= args.tol
                report["pass"] = report["pass"] and ok
                report["migration"].append(
                    {"i": i, "j": j, "ba3_mean": mean, "ba3_sd": sd,
                     "realized": realized, "z": z, "within_tol": ok})
                print("  %-8s %-14s %-14.4f %-8.2f %s"
                      % ("%d<-%d" % (i, j), "%.4f(%.4f)" % (mean, sd), realized, z,
                         "yes" if ok else "NO"))
            else:
                print("  %-8s %-14s %-14.4f  (no estimate)" % ("%d<-%d" % (i, j), "-", realized))

    # ---- dispersal sex bias: phiMove (movement), phiBreed (gene flow), rho ----
    prm = truth["params"]
    smry = truth["summary"]
    targets = [
        ("phiMove",  prm.get("phi_move"),  smry.get("empirical_phi_move"),  "movement (age-1)"),
        ("phiBreed", prm.get("phi_breed"), smry.get("empirical_phi_breed"), "gene flow (age-2)"),
        ("rho",      None,                 smry.get("empirical_rho"),       "residents"),
    ]
    print("\nSex-biased dispersal  (female fractions):")
    print("  %-9s %-16s %-10s %-8s %s" % ("param", "BA3 mean(SD)", "target", "z", "within tol"))
    for key, true_val, emp_val, label in targets:
        target = true_val if true_val is not None else emp_val  # rho has no input param
        est = sexbias.get(key)
        if est is None:
            print("  %-9s (not found in output)   [%s]" % (key, label))
            continue
        mean, sd = est
        z = (target - mean) / sd if (target is not None and sd > 0) else float("nan")
        ok = (target is None) or abs(z) <= args.tol
        if true_val is not None:
            report["pass"] = report["pass"] and ok
        report["sexbias"][key] = {"ba3_mean": mean, "ba3_sd": sd, "target": target,
                                  "z": z, "within_tol": ok, "label": label}
        tstr = ("%.4f" % target) if target is not None else "-"
        print("  %-9s %-16s %-10s %-8.2f %s   [%s]"
              % (key, "%.4f(%.4f)" % (mean, sd), tstr, z,
                 "yes" if ok else "NO", label))

    print("\n" + "-" * 64)
    print("OVERALL: %s" % ("PASS" if report["pass"] else "FAIL"))

    if args.json:
        print(json.dumps(report, indent=1))

    sys.exit(0 if report["pass"] else 1)


if __name__ == "__main__":
    main()
