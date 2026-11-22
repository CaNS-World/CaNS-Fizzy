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
    for field in ["u","v","w","p","psi","s_001","s_002","s_003"]:
        assert "fld_"+field+".bin" in names, f"Checkpoint of {field} not found."
    for name in names:
        data_one_step = np.fromfile("data-one-step/"+name)
        data_two_step = np.fromfile("data-two-step/"+name)
        assert np.all(np.isfinite(data_one_step)), f"{name}: invalid values."
        assert np.array_equal(data_one_step,data_two_step), \
               f"{name}: restarted run differs, by at most {np.max(np.abs(data_one_step-data_two_step))}."
    #
    # three scalars are transported: the second one has the inputs of the first, and must be equal to it; the
    # third one has other properties, boundary values and a source, and must differ
    #
    s1,s2,s3 = (np.fromfile(f"data-one-step/fld_s_00{i}.bin")[:-2] for i in (1,2,3))
    assert np.array_equal(s1,s2), "Two scalars with the same inputs differ."
    assert np.max(np.abs(s3-s1)) > 1.e-3, "The third scalar does not differ from the first one."
if __name__ == "__main__":
    test_restart_two_fluid()
    print("Passed!")
