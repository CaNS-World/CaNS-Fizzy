#!/usr/bin/env python
def test_marangoni_markers():
    import numpy as np
    #
    # thermocapillary migration of drops in 2D (the example `thermocapillary_migration_2d` on 32^2 cells, 100 steps),
    # with the geometric method:
    #   one       : a single marker, one drop
    #   second    : two markers, the drop in the second one; must reproduce the case `one` exactly
    #   pair-same : a single marker, two smaller drops
    #   pair-split: two markers, one drop each; must reproduce the case `pair-same` exactly, the drops being far apart
    #
    def load(case,field):
        data = np.fromfile("data-"+case+"/fld_"+field+".bin")
        assert np.all(np.isfinite(data)), f"{case}, {field}: invalid values."
        return data[:-2] # the last two entries are the time and the time step number
    for ref,case in [("one","second"),("pair-same","pair-split")]:
        for field in ["u","v","w","p","psi","s_001"]:
            assert np.array_equal(load(ref,field),load(case,field)), \
                   f"{field}: the case `{case}` differs from the case `{ref}`."
        #
        # the drops must be moving (towards the hot wall, along z)
        #
        psi = load(ref,"psi"); w = load(ref,"w")
        w_drop = np.sum(psi*w)/np.sum(psi)
        assert w_drop > 1.e-3, f"The drops of the case `{ref}` are not migrating: mean velocity {w_drop}."
if __name__ == "__main__":
    test_marangoni_markers()
    print("Passed!")
