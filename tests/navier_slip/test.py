#!/usr/bin/env python
def test_navier_slip():
    import numpy as np
    #
    # Navier slip at a wall: plane channel of height 1 driven by a unit body force, with a slip length of 0.1
    # at the lower wall and no slip at the upper one, single fluid with unit viscosity; at steady state,
    #   u(z) = (-z^2 + a*z + lambda*a)/2, with a = 1/(1+lambda)
    #
    nz, slip = 32, 0.1
    u = np.fromfile('data/vey_fld_0040000.bin').reshape((8, nz), order='F')
    z = (np.arange(nz)+0.5)/nz
    a = 1./(1.+slip)
    u_exact = 0.5*(-z**2+a*z+slip*a)
    err = np.max(np.abs(u-u_exact[None,:]))/np.max(u_exact)
    assert err < 2.e-3, f"Velocity profile differs from the exact one by {err} (about 8.e-4 is expected)."
if __name__ == "__main__":
    test_navier_slip()
    print("Passed!")
