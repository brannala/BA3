#!/usr/bin/env python3
"""Integrated autocorrelation time (IAT) and effective sample size (ESS) for
BA3 trace columns. ESS = N / tau, tau = 1 + 2*sum_{t>=1} rho(t), summed with the
initial-positive-sequence (Geyer) truncation. Reports per-migration-column ESS."""
import sys, numpy as np

def iat(x):
    x = np.asarray(x, float); x = x - x.mean()
    n = len(x); v = np.dot(x, x) / n
    if v == 0: return 1.0
    # autocovariance via FFT
    f = np.fft.rfft(x, n=2*n)
    acov = np.fft.irfft(f*np.conj(f))[:n].real / n
    rho = acov / acov[0]
    # Geyer initial positive sequence: sum pairs until a pair sum <= 0
    tau = 1.0; t = 1
    while t + 1 < n:
        pair = rho[t] + rho[t+1]
        if pair <= 0: break
        tau += 2.0 * pair; t += 2
    return max(tau, 1.0)

def main():
    label, fn, burn = sys.argv[1], sys.argv[2], int(sys.argv[3])
    hdr = open(fn).readline().rstrip('\n').split('\t')
    mig_cols = [i for i,h in enumerate(hdr) if h.startswith('m[')]
    data = np.loadtxt(fn, skiprows=1+burn, usecols=mig_cols)
    if data.ndim == 1: data = data[:,None]
    n = data.shape[0]
    taus, esss = [], []
    for j in range(data.shape[1]):
        col = data[:,j]
        if col.std() == 0: continue
        t = iat(col); taus.append(t); esss.append(n/t)
    taus, esss = np.array(taus), np.array(esss)
    print("%-9s n=%d  cols=%d  IAT median=%.1f (min %.1f max %.1f)  ESS median=%.0f (min %.0f)  ESS/iter=%.4f"
          % (label, n, len(taus), np.median(taus), taus.min(), taus.max(),
             np.median(esss), esss.min(), np.median(esss)/n))

if __name__ == "__main__": main()
