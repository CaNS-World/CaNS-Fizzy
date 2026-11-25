#!/usr/bin/env python
def test_options():
    import numpy as np
    #
    # rising bubble of the test `rising_bubble` after 1000 time steps (t = 0.5, when its centre of mass has
    # risen from 0.5 to about 0.553), with options that no other test selects, each against its reference case:
    #   acdi-nosdf, thinc-nosdf: normals and curvature from the phase indicator (`is_sdf_normals = F`)
    #   acdi-nobrack           : surface tension force without the density scaling
    #   acdi-harmonic          : harmonic interpolation of the viscosity
    #   plic-hfnormals         : height-function normals, with the transport along with the sweeps
    #   plic-nolimiter         : geometric method without the viscosity limiter
    # the centre of mass must agree with that of the reference case within a few times the difference
    # measured when the test was written (given below), and the fields must differ: the option acts;
    # note that ACDI without the signed-distance function is the one that is far from its reference
    #
    ny, nz = 64, 128
    z = (np.arange(nz)+0.5)*2./nz
    def load(case,field):
        data = np.fromfile("data-"+case+"/fld_"+field+".bin")
        assert np.all(np.isfinite(data)), f"{case}, {field}: invalid values."
        return data[:-2].reshape((ny,nz),order="F") # the last two entries are the time and the time step number
    def centre(case):
        gas = 1.-load(case,"psi")
        return np.sum(gas*z[None,:])/np.sum(gas)
    for case, ref, tol in [("acdi-nosdf"    ,"acdi" ,8.e-3),  # 3.7e-3
                           ("acdi-nobrack"  ,"acdi" ,2.e-4),  # 2.3e-5
                           ("acdi-harmonic" ,"acdi" ,3.e-3),  # 9.0e-4
                           ("thinc-nosdf"   ,"thinc",2.e-4),  # 1.6e-5
                           ("plic-hfnormals","plic" ,2.e-4),  # 2.0e-6
                           ("plic-nolimiter","plic" ,1.e-3)]: # 2.5e-4
        diff = abs(centre(case)-centre(ref))
        assert diff < tol, f"{case}: centre of mass {centre(case)}, against {centre(ref)} of the case {ref}."
        assert np.max(np.abs(load(case,"psi")-load(ref,"psi"))) > 1.e-6, f"{case}: same fields as the case {ref}."
    for case in ["acdi","thinc","plic"]:
        assert abs(centre(case)-0.5535) < 2.e-3, f"{case}: centre of mass {centre(case)}, expected about 0.5535."
if __name__ == "__main__":
    test_options()
    print("Passed!")
