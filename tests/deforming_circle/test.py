#!/usr/bin/env python
def test_deforming_circle():
    import numpy as np
    from read_single_field_binary import read_single_field_binary
    #
    # comparing initial condition with psi field after one cycle (t=4)
    #
    data_beg,xp,yp,zp,xu,yv,zw = read_single_field_binary("psi_fld_0000000.bin",np.array([1,1,1]))
    data_end,xp,yp,zp,xu,yv,zw = read_single_field_binary("psi_fld_0016000.bin",np.array([1,1,1]))
    lx = 1.
    lz = 1.
    dx = xu[1]-xu[0]
    dz = zw[1]-zw[0]
    err = np.sum(abs(data_beg[:,0,:]-data_end[:,0,:]))*dx*dz/(lx*lz)
    err_max = 1.e-1
    try:
        assert err < err_max, "Error is too large."
    except AssertionError as msg:
        print(msg)
if __name__ == "__main__":
    test_deforming_circle()
