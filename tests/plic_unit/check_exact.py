#!/usr/bin/env python3
#
# checks the plane volumes written by test_plic (volume_samples.txt) against exact rational arithmetic
#
import sys
from fractions import Fraction
from itertools import product

def vol_exact(m, alpha):
    a = Fraction(alpha) - sum(min(Fraction(q), 0) for q in m)
    b = [abs(Fraction(q)) for q in m if q != 0.]
    if not b:
        return Fraction(int(a > 0))
    v, den, p = Fraction(0), 1, len(b)
    for q in b:
        den *= q
    for s in product((0, 1), repeat=p):
        t = a - sum(si*bi for si, bi in zip(s, b))
        if t > 0:
            v += (-1)**sum(s)*t**p
    fact = {1: 1, 2: 2, 3: 6}[p]
    return min(max(v/(fact*den), Fraction(0)), Fraction(1))

errmax = 0.
for line in open(sys.argv[1] if len(sys.argv) > 1 else 'volume_samples.txt'):
    mx, my, mz, alpha, vol = (float(x) for x in line.split())
    errmax = max(errmax, abs(float(Fraction(vol) - vol_exact((mx, my, mz), alpha))))
print(f' volume: max abs error vs exact rational arithmetic (all normal kinds): {errmax:.2e}')
if errmax > 1e-14:
    print('  FAILED: volume error above 1e-14')
    sys.exit(1)
