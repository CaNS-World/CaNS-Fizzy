#!/usr/bin/env python
def test_single_precision():
    import numpy as np
    #
    # single- against double-precision runs of a few short cases, one for each interface method and for the
    # features that depend most on round-off (height functions, contact angle, several markers, Marangoni
    # stresses): the largest difference of the final velocity, pressure (without its mean) and phase indicator
    # must stay below the bounds, which are about ten times the differences measured when the test was written
    # (absolute differences: the test prints them, along with the largest velocity and pressure of each case)
    #
    #                     velocity  pressure  phase indicator
    bounds = {"bubble-acdi"      : (5.e-4, 2.e+0 ,1.e-3),
              "bubble-thinc"     : (5.e-4, 2.e+0 ,1.e-3),
              "bubble-plic"      : (6.e-3, 2.e+1 ,6.e-3),
              "bubble-thincsplit": (2.e-3, 1.e+0 ,2.e-3),
              "drop-plic"        : (5.e-5, 2.e-4 ,1.e-3),
              "advection-plic"   : (1.e-5, 1.e-12,3.e-1), # (no flow solution; thin filaments: large local differences)
              "markers-plic"     : (1.e-4, 3.e-3 ,6.e-4),
              "marangoni-plic"   : (3.e-4, 7.e-3 ,2.e-4)}
    def load(case,precision,field):
        kind = {"double": np.float64, "single": np.float32}[precision]
        data = np.fromfile("data-"+case+"-"+precision+"/fld_"+field+".bin",dtype=kind)
        assert np.all(np.isfinite(data)), f"{case}, {precision} precision, {field}: invalid values."
        return data[:-2].astype(np.float64) # the last two entries are the time and the time step number
    def difference(case,field):
        a = load(case,"double",field)
        b = load(case,"single",field)
        if field == "p": # (the pressure is defined up to a constant)
            a -= np.mean(a)
            b -= np.mean(b)
        return np.max(np.abs(a-b))
    failed = []
    for case,(vel_max,pre_max,psi_max) in bounds.items():
        vel = max(difference(case,field) for field in ["u","v","w"])
        pre = difference(case,"p")
        psi = difference(case,"psi")
        vel_ref = max(np.max(np.abs(load(case,"double",field))) for field in ["u","v","w"])
        pre_ref = np.max(np.abs(load(case,"double","p")-np.mean(load(case,"double","p"))))
        print(f"{case:18s}: differences in velocity, pressure and phase indicator: {vel:.1e} {pre:.1e} {psi:.1e}"
              f" (largest velocity and pressure: {vel_ref:.1e} {pre_ref:.1e})")
        if not (vel <= vel_max and pre <= pre_max and psi <= psi_max):
            failed.append(f"{case}: {vel:.1e} {pre:.1e} {psi:.1e} against the bounds {vel_max:.1e} {pre_max:.1e} {psi_max:.1e}")
    assert not failed, "single and double precision differ by more than the bounds: "+"; ".join(failed)
if __name__ == "__main__":
    test_single_precision()
    print("Passed!")
