#!/usr/bin/env python
def test_single_fluid():
    import numpy as np
    #
    # lid-driven cavity of a single fluid, 200 time steps: the velocity and the pressure of the run with
    # phase 1 everywhere must be those of the run with phase 2 everywhere, to round-off, although the
    # properties of the other phase differ; the pressure is compared without its mean
    #
    for transport in ['ab2', 'one-step', 'sweeps']:
        for field in ['u', 'v', 'w', 'p']:
            a = np.fromfile(f'data-{transport}-1/fld_{field}.bin')[:-2]
            b = np.fromfile(f'data-{transport}-2/fld_{field}.bin')[:-2]
            assert np.all(np.isfinite(a)), f"{transport}, {field}: invalid values."
            if field == 'p':
                a, b = a-np.mean(a), b-np.mean(b)
            err = np.max(np.abs(a-b))
            assert err < 1.e-11, f"{transport}, {field}: phase 1 and phase 2 differ by {err}."
        v = np.fromfile(f'data-{transport}-1/fld_v.bin')[:-2]
        assert np.max(np.abs(v)) > 0.1, f"{transport}: the flow did not develop."
if __name__ == "__main__":
    test_single_fluid()
    print("Passed!")
