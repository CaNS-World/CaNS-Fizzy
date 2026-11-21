#!/usr/bin/env python
def test_plic_advection():
    import numpy as np
    #
    # geometric (PLIC) volume-of-fluid method with prescribed, time-reversed velocity fields, at a CFL number of 0.25:
    #   2d: circle in the serpentine vortex (period 4), 64^2 cells
    #   3d: sphere in the deformation field of LeVeque (period 3), 32^3 cells
    # the interface must return to its initial shape, with the volume conserved and the volume fraction bounded
    # to round-off; the tolerances of the shape error (L1 norm) are about 1.3 times the errors of the method
    #
    cases = {"2d": (1024,3.0e-3), "3d": (768,1.0e-2)}
    for case,(nstep,err_max) in cases.items():
        psi_beg = np.fromfile("data-"+case+"/psi_fld_0000000.bin")
        psi_end = np.fromfile("data-"+case+f"/psi_fld_{nstep:07d}.bin")
        assert np.all(np.isfinite(psi_end)), f"{case}: invalid values."
        err = np.mean(np.abs(psi_end-psi_beg)) # the cells have the same size
        assert err < err_max, f"{case}: the shape error is too large: {err} >= {err_max}."
        #
        # log file: time step, time, volume, relative volume change, clipped volume,
        #           extrema of the volume fraction before clipping, CFL numbers, number of sub-steps
        #
        log = np.loadtxt("data-"+case+"/log_vof.out")
        dvol = np.max(np.abs(log[:,3]))
        assert dvol < 1.e-12, f"{case}: the volume is not conserved: {dvol}."
        assert np.min(log[:,5]) > -1.e-12 and np.max(log[:,6]) < 1.+1.e-12, f"{case}: the volume fraction is not bounded."
if __name__ == "__main__":
    test_plic_advection()
    print("Passed!")
