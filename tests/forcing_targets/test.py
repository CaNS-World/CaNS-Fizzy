#!/usr/bin/env python
def test_forcing_targets():
    import os
    import numpy as np
    case = os.environ.get('FORCING_CASE', 'conduction')
    if case == 'conduction':
        #
        # a fluid at rest between two walls at s = 0, with unit rho*cp and conductivity ka = 0.1, starts at s = 1
        # (`scalf`) and is held at that mean by the forcing: at steady state s = 6*z*(1-z), and the uniform source
        # is 12*ka; the discrete problem with nz cells of size h has the exact solution s = q/(2*ka)*(z*(1-z)+h**2/4),
        # with the source q = 12*ka/(1+2*h**2)
        #
        nz, ka = 16, 0.1
        s = np.fromfile('data/sca_001_fld_0001000.bin').reshape((8, nz), order='F')
        z = (np.arange(nz)+0.5)/nz
        q = 12.*ka/(1.+2./nz**2)
        err = np.max(np.abs(s-q/(2.*ka)*(z*(1.-z)+0.25/nz**2)[None, :]))
        assert err < 1.e-6, f"Scalar profile differs from the exact one by {err}."
        log = np.loadtxt('data/forcing_scalar_001.out')
        source, mean = log[-1, 1], log[:, 4]
        err = abs(source/q-1.)
        assert err < 1.e-6, f"Source of the forcing differs from the steady one by {err} (relative)."
        err = np.max(np.abs(mean-1.))
        assert err < 1.e-12, f"Mean of the scalar departs from its target by {err}."
    elif case == 'bubble':
        #
        # a rising bubble with rho*cp = 400 in a fluid with rho*cp = 1000, conductivities 0.05 and 0.5, between a
        # wall at s = 0 and one at s = 2, from s = 0.5 (`scalf`), on a grid stretched along z: the mean of rho*cp*s
        # must stay constant, i.e., its measured rate of change must vanish next to the source
        #
        log = np.loadtxt('data/forcing_scalar_001.out')
        source, residual, total, mean = log[:, 1], log[:, 2], log[:, 3], log[:, 4]
        err = np.max(np.abs(residual))/np.max(np.abs(source))
        assert err < 1.e-6, f"Residual of the balance of the scalar is {err} of the source."
        err = np.max(np.abs(total/total[0]-1.))
        assert err < 1.e-10, f"Mean of rho*cp*s changes by {err} (relative)."
        assert abs(mean[0]-0.5) < 1.e-6, f"Initial weighted mean of the scalar is {mean[0]}, not the target 0.5."
    else:
        #
        # a bubble ten times lighter than the fluid around it in a forced channel flow: the mass-averaged
        # velocity (total momentum over total mass) of the initial condition is `velf` = 0.7, and it is kept
        #
        n, velf = (12, 12, 12), 0.7
        dzf = np.loadtxt('data/grid.out')[1:n[2]+1, 3]
        for step, tol in (('0000000', 1.e-12), ('0000040', 1.e-6)):
            u, psi = (np.fromfile(f'data/{q}_fld_{step}.bin').reshape(n, order='F') for q in ('vex', 'psi'))
            rho = 1.+9.*0.5*(psi+np.roll(psi, -1, axis=0))
            vel = np.sum(rho*u*dzf[None, None, :])/np.sum(rho*dzf[None, None, :])
            assert abs(vel/velf-1.) < tol, f"Mass-averaged velocity at step {int(step)} is {vel}, not the target {velf}."
if __name__ == "__main__":
    test_forcing_targets()
    print("Passed!")
