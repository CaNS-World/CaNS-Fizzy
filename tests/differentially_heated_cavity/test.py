#!/usr/bin/env python
def test_dht():
    import numpy as np
    from read_single_field_binary import read_single_field_binary
    #
    # testing u
    #
    data,xp,yp,zp,xu,yv,zw = read_single_field_binary("sca_fld_0050000.bin",np.array([1,1,1]))
    islice = np.size(data[0,0,:])//2
    tw = -0.5
    l  = 1.
    nusselt = (1./l)*np.sum((data[0,0,:]-tw)/(2.*xp[0])*(2.*zp[0]))*l/(-tw)
    data,xp,yp,zp,xu,yv,zw = read_single_field_binary("vex_fld_0050000.bin",np.array([1,1,1]))
    islice = int(np.size(data[0,0,:])/2)
    umax = np.max(data[islice,0,:]+data[islice+1,0,:])/2
    data,xp,yp,zp,xu,yv,zw = read_single_field_binary("vez_fld_0050000.bin",np.array([1,1,1]))
    islice = int(np.size(data[0,0,:])/2)
    wmax = np.max(data[:,0,islice]+data[:,0,islice+1])/2
    #
    nusselt_ref = 8.8252
    np.testing.assert_allclose([nusselt], [nusselt_ref], rtol=1.e-2, atol=1.e-2)
if __name__ == "__main__":
    test_dht()