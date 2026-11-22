#!/usr/bin/env python
def test_implicit_diffusion():
    import numpy as np
    #
    # a viscous two-fluid flow with a scalar, on a grid stretched along z:
    # the runs with implicit diffusion, along z (`is_impdiff_1d = T`) and in all directions with constant
    # coefficients (`is_impdiff = T`), at time steps 10 to 20 times larger than the limit of the explicit
    # scheme, must agree with an explicit run at a small time step;
    # the tolerances are a few times the differences that the second-order schemes yield here; those of the
    # geometric volume-of-fluid method (variable TEST_TWO_FLUID_INPUT, see `tests/testall.sh`) are larger with
    # constant coefficients, whose splitting error is larger across a sharp interface
    #
    import os
    is_sharp = "plic" in os.environ.get("TEST_TWO_FLUID_INPUT","")
    def load(case,field):
        data = np.fromfile("data-"+case+"/fld_"+field+".bin")
        assert np.all(np.isfinite(data)), f"{case}, {field}: invalid values."
        return data[:-2] # the last two entries are the time and the time step number
    vel_scale = max(np.max(np.abs(load("explicit",field))) for field in ["u","v","w"])
    assert vel_scale > 0.
    tols = {"u": 2.e-3*vel_scale, "v": 2.e-3*vel_scale, "w": 2.e-3*vel_scale, "psi": 1.e-5, "s": 1.e-4}
    for case in ["implicit-z","implicit-3d"]:
        for field,tol in tols.items():
            if is_sharp and case == "implicit-3d": tol = 2.5*tol
            err = np.max(np.abs(load(case,field)-load("explicit",field)))
            assert err < tol, f"{case}, {field}: implicit and explicit runs differ by {err} (tolerance: {tol})."
if __name__ == "__main__":
    test_implicit_diffusion()
    print("Passed!")
