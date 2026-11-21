#!/usr/bin/env python
def test_plic_momentum_drop():
    import numpy as np
    #
    # geometric volume-of-fluid method with the momentum and the scalar transported along with the sweeps
    # of the volume fraction (`is_vof_sweep_transport = T`): inviscid drop in a uniform diagonal flow, with
    # density and heat capacity ratios of 1e6, on 64^2 cells at a CFL number of 0.25; the velocity and the
    # uniform scalar must stay uniform to round-off (mass, momentum and energy are transported consistently),
    # and the volume fraction conserved and bounded
    #
    ng = np.array([64, 2, 64])
    u = np.fromfile('data-drop/vex_fld_0000256.bin').reshape(ng, order='F')
    w = np.fromfile('data-drop/vez_fld_0000256.bin').reshape(ng, order='F')
    s = np.fromfile('data-drop/sca_fld_0000256.bin').reshape(ng, order='F')
    err_max = 1.e-9
    err = max(np.max(abs(u-1.)), np.max(abs(w-1.)))
    assert err < err_max, f"Uniform flow not preserved: {err} >= {err_max}."
    err = np.max(abs(s-1.))
    assert err < err_max, f"Uniform scalar not preserved: {err} >= {err_max}."
    log = np.loadtxt('data-drop/log_vof.out')
    assert np.max(abs(log[:,3])) < 1.e-12, f"Volume not conserved: {np.max(abs(log[:,3]))}."
    assert np.min(log[:,5]) > -1.e-12 and np.max(log[:,6]) < 1.+1.e-12, "Volume fraction out of bounds."
def test_plic_momentum_restart():
    import os
    import glob
    import numpy as np
    #
    # rising bubble with a scalar (the case of `restart_two_fluid`) and the same transport: a run
    # restarted halfway with `is_full_restart = T` must reproduce the uninterrupted run exactly
    #
    names = sorted(os.path.basename(f) for f in glob.glob("data-restart-one-step/fld_*.bin"))
    for field in ["u","v","w","p","psi","s","uadv","vadv","wadv","sadv"]:
        assert "fld_"+field+".bin" in names, f"Checkpoint of {field} not found."
    for name in names:
        data_one_step = np.fromfile("data-restart-one-step/"+name)
        data_two_step = np.fromfile("data-restart-two-step/"+name)
        assert np.all(np.isfinite(data_one_step)), f"{name}: invalid values."
        assert np.array_equal(data_one_step,data_two_step), \
               f"{name}: restarted run differs, by at most {np.max(np.abs(data_one_step-data_two_step))}."
if __name__ == "__main__":
    test_plic_momentum_drop()
    test_plic_momentum_restart()
    print("Passed!")
