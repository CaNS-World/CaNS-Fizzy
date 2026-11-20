#!/usr/bin/env python
def test_restart_two_fluid():
    import os
    import glob
    import numpy as np
    #
    # a run restarted halfway with `is_full_restart = T` must reproduce the uninterrupted run exactly:
    # all checkpoint fields are compared bitwise
    #
    names = sorted(os.path.basename(f) for f in glob.glob("data-one-step/fld_*.bin"))
    for field in ["u","v","w","p","psi","s"]:
        assert "fld_"+field+".bin" in names, f"Checkpoint of {field} not found."
    for name in names:
        data_one_step = np.fromfile("data-one-step/"+name)
        data_two_step = np.fromfile("data-two-step/"+name)
        assert np.all(np.isfinite(data_one_step)), f"{name}: invalid values."
        assert np.array_equal(data_one_step,data_two_step), \
               f"{name}: restarted run differs, by at most {np.max(np.abs(data_one_step-data_two_step))}."
if __name__ == "__main__":
    test_restart_two_fluid()
    print("Passed!")
