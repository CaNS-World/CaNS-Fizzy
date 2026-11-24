#!/usr/bin/env python
def test_markers():
    import numpy as np
    #
    # two drops at rest, half a grid cell apart, with surface tension (64^2 cells):
    #   one  : a single marker; the drops merge, which sets the fluid in motion
    #   same : two markers, both drops in the first one; must reproduce the case `one` exactly
    #   split: one drop per marker; the drops must stay apart and at rest, each marker conserving its volume
    #
    n = 64
    def load(case,field):
        data = np.fromfile("data-"+case+"/fld_"+field+".bin")
        assert np.all(np.isfinite(data)), f"{case}, {field}: invalid values."
        return data[:-2] # the last two entries are the time and the time step number
    for field in ["u","v","w","p","psi"]:
        assert np.array_equal(load("one",field),load("same",field)), \
               f"{field}: two markers with all the drops in the first one differ from a single marker."
    #
    # largest velocity, and width of the bridge between the two drops (in grid cells)
    #
    def vel_max(case):
        return max(np.max(np.abs(load(case,field))) for field in ["v","w"])
    def bridge(case):
        return np.mean(np.sum(load(case,"psi").reshape((n,n),order="F")[n//2-1:n//2+1,:],axis=1))
    assert vel_max("one"  ) > 1.e-1, f"The drops of a single marker did not merge: largest velocity {vel_max('one')}."
    assert vel_max("split") < 1.e-2, f"The drops of two markers did not stay at rest: largest velocity {vel_max('split')}."
    assert bridge("one") > 1.5*bridge("split"), f"Bridge between the drops: {bridge('one')} against {bridge('split')} cells."
    #
    # volume of each marker, against that of the initial condition
    #
    for marker in ["001","002"]:
        vol_beg = np.sum(np.fromfile("data-split/psi_"+marker+"_fld_0000000.bin"))
        vol_end = np.sum(np.fromfile("data-split/psi_"+marker+"_fld_0000600.bin"))
        assert abs(vol_end/vol_beg-1.) < 1.e-10, f"Marker {marker} does not conserve its volume: {vol_end/vol_beg-1.}."
if __name__ == "__main__":
    test_markers()
    print("Passed!")
