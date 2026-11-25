#!/usr/bin/env python
def test_forced_hit():
    import os
    import numpy as np
    #
    # large-scale forcing (`is_forced_hit = T`) of a single fluid at rest in a periodic box of side 2*pi, with
    # viscosity nu; with k = 2 and f0 = `forced_hit_amp` = 0.5 (see `src/forcing.f90`):
    #
    #   'abc': f is a Beltrami field (curl f = k f), so u = f*(1-exp(-nu*k^2*t))/(nu*k^2) is an exact solution
    #   'tgv': u = f*(1-exp(-3*nu*k^2*t))/(3*nu*k^2), plus the effect of advection, of relative size f0*k*t^2
    #
    ftype = os.environ.get('FORCED_HIT_TYPE', 'tgv')
    n, nu, k, f0, time = (32, 24, 24), 0.05, 2., 0.5, 0.1
    u, v, w = (np.fromfile(f'data/ve{q}_fld_0000100.bin').reshape(n, order='F') for q in 'xyz')
    xc, yc, zc = ((np.arange(m)+0.5)*2.*np.pi/m for m in n)
    xf, yf, zf = ((np.arange(m)+1.0)*2.*np.pi/m for m in n)
    def grid(x, y, z): return np.meshgrid(x, y, z, indexing='ij')
    if ftype == 'abc':
        gain = f0*(1.-np.exp(-nu*k**2*time))/(nu*k**2)
        x, y, z = grid(xf, yc, zc); u_ref = gain*(np.sin(k*z)+np.cos(k*y))
        x, y, z = grid(xc, yf, zc); v_ref = gain*(np.sin(k*x)+np.cos(k*z))
        x, y, z = grid(xc, yc, zf); w_ref = gain*(np.sin(k*y)+np.cos(k*x))
        tol, expected = 1.e-3, '4.e-4'
    else:
        gain = f0*(1.-np.exp(-3.*nu*k**2*time))/(3.*nu*k**2)
        x, y, z = grid(xf, yc, zc); u_ref =  gain*np.sin(k*x)*np.cos(k*y)*np.cos(k*z)
        x, y, z = grid(xc, yf, zc); v_ref = -gain*np.cos(k*x)*np.sin(k*y)*np.cos(k*z)
        w_ref = np.zeros(n)
        tol, expected = 5.e-3, '2.e-3'
    err = max(np.max(np.abs(q-q_ref)) for q, q_ref in ((u, u_ref), (v, v_ref), (w, w_ref)))/gain
    assert err < tol, f"Velocity differs from the reference by {err} of its amplitude (about {expected} is expected)."
if __name__ == "__main__":
    test_forced_hit()
    print("Passed!")
