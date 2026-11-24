#!/usr/bin/env python
def test_markers():
    import numpy as np
    #
    # two drops at rest, two grid cells apart, with surface tension (64^2 cells, ACDI):
    #   one  : a single marker; the drops merge
    #   same : two markers, both drops in the first one; must reproduce the case `one` exactly
    #   split: one drop per marker; the drops must stay apart, each marker conserving its volume
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
    # phase indicator between the two drops
    #
    def psi_mid(case):
        return np.mean(load(case,"psi").reshape((n,n),order="F")[n//2-1:n//2+1,n//2-1:n//2+1])
    assert psi_mid("one"  ) > 0.9, f"The drops of a single marker did not merge: {psi_mid('one')}."
    assert psi_mid("split") < 0.6, f"The drops of two markers merged: {psi_mid('split')}."
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
