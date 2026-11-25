#!/usr/bin/env python
def test_deforming_circle():
    import numpy as np
    from read_single_field_binary import read_single_field_binary
    #
    # comparing initial condition with psi field after one cycle (t=4): L1 norm of the difference, over the
    # area of the circle (0.010 with ACDI and 0.005 with THINC/QQ on this grid of 256 x 256 cells)
    #
    data_beg,xp,yp,zp,xu,yv,zw = read_single_field_binary("data/","psi_fld_0000000.bin",np.array([1,1,1]))
    data_end,xp,yp,zp,xu,yv,zw = read_single_field_binary("data/","psi_fld_0016000.bin",np.array([1,1,1]))
    dy = yv[1]-yv[0]
    dz = zw[1]-zw[0]
    area = np.sum(1.-data_beg[0,:,:])*dy*dz
    err = np.sum(abs(data_beg[0,:,:]-data_end[0,:,:]))*dy*dz/area
    err_max = 3.e-2
    assert err < err_max, f"Error is too large: {err} >= {err_max}."
if __name__ == "__main__":
    test_deforming_circle()
