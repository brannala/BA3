#!/usr/bin/env python3
"""Phase-3 posterior-equivalence gate for the collapsed sampler.

Runs BA3 in standard and --collapse modes on the SAME base-model dataset (2-column
metadata, no sex) and the SAME seed, and checks that the migration-rate posteriors
agree within MC error (max |m_std - m_col| / combined SD below a threshold).

Usage:
  python3 sim/check_collapsed_posterior.py ./BA3 data.vcf data_meta2.txt [iters]
"""
import sys, re, subprocess, tempfile, os

THRESH = 0.6   # max |z| considered "agree" (SD-scaled)


def run(ba3, vcf, meta, iters, collapse):
    out = tempfile.mktemp(suffix=".txt")
    cmd = [ba3, '-V', vcf, '-M', meta, '-i', str(iters), '-b', str(iters // 5),
           '-n', '200', '-s', '17', '-o', out]
    if collapse:
        cmd.insert(1, '-c')
    subprocess.run(cmd, capture_output=True, text=True)
    txt = open(out).read()
    os.remove(out)
    return txt


def matrix(txt):
    L = txt.splitlines()
    out = {}
    for k, l in enumerate(L):
        if 'Migration Rate Matrix' in l:
            r = k + 3
            while r < len(L) and re.match(r'\s*\[\d+\]', L[r]):
                i = int(re.search(r'\[(\d+)\]', L[r]).group(1))
                for j, (m, s) in enumerate(re.findall(r'([0-9.]+)\(([0-9.]+)\)', L[r])):
                    out[(i, j)] = (float(m), float(s))
                r += 1
            break
    return out


def main():
    ba3, vcf, meta = sys.argv[1], sys.argv[2], sys.argv[3]
    iters = int(sys.argv[4]) if len(sys.argv) > 4 else 500000
    a = matrix(run(ba3, vcf, meta, iters, False))
    b = matrix(run(ba3, vcf, meta, iters, True))
    n = max(i for i, _ in a) + 1
    mx = 0.0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            (ma, sa), (mb, sb) = a[(i, j)], b[(i, j)]
            z = abs(ma - mb) / ((sa + sb) / 2 if (sa + sb) else 1)
            mx = max(mx, z)
            print("  m[%d<-%d]  std %.4f(%.4f)  collapsed %.4f(%.4f)  z=%.2f"
                  % (i, j, ma, sa, mb, sb, z))
    ok = mx < THRESH
    print("max |z| = %.2f  ->  %s" % (mx, "PASS (posteriors agree)" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
