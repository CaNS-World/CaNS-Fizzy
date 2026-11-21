#!/usr/bin/env python
def test_viscosity_limiter():
    import numpy as np
    #
    # start of the rising bubble of Hysing et al., case 2 (density ratio 1000, viscosity ratio 100), with the
    # sharp interface of the geometric volume-of-fluid method: the viscosity limiter (`is_visc_limiter`, on
    # by default with this method) keeps the explicit diffusion stable next to the interface; without it,
    # the velocity is 1000 times larger after these 60 time steps
    #
    vel_max = max(np.max(np.abs(np.fromfile('data/ve'+c+'_fld_0000060.bin'))) for c in ['x','y','z'])
    assert np.isfinite(vel_max), "Invalid velocity."
    assert vel_max < 0.05, f"Velocity too large: {vel_max} (about 0.013 is expected)."
if __name__ == "__main__":
    test_viscosity_limiter()
    print("Passed!")
