#!/usr/bin/env python
def test_ldc():
    import numpy as np
    from read_single_field_binary import read_single_field_binary
    #
    # velocity profiles along the two centrelines of the cavity at the end of the run, against the values
    # tabulated by Ghia et al., J. Comput. Phys. 48 (1982) 387-411, for a Reynolds number of 100;
    # each velocity component is taken at the cell faces that lie on the centreline (y = 0.5 or z = 0.5),
    # and the values at the walls are added for the interpolation to the tabulated points;
    # the largest differences are 0.004 and 0.009 on this grid of 64 x 64 cells
    #
    err_max = 2.e-2
    #
    # velocity along y (that of the lid), as a function of z
    #
    data,xp,yp,zp,xu,yv,zw = read_single_field_binary("data/","vey_fld_0020000.bin",np.array([1,1,1]))
    islice = np.size(data[0,:,0])//2-1
    assert abs(yv[islice]-0.5) < 1.e-12, f"Faces of the velocity along y not on the centreline: {yv[islice]}."
    data_ref = np.loadtxt("ghiau.txt") # https://gist.github.com/ivan-pi/3e9326d18a366ffe6a8e5bfda6353219
    data_interp = np.interp(data_ref[:,0],np.r_[0.,zp,1.],np.r_[0.,data[0,islice,:],1.])
    np.testing.assert_allclose(data_interp, data_ref[:,1], rtol=0., atol=err_max)
    #
    # velocity along z, as a function of y
    #
    data,xp,yp,zp,xu,yv,zw = read_single_field_binary("data/","vez_fld_0020000.bin",np.array([1,1,1]))
    islice = np.size(data[0,0,:])//2-1
    assert abs(zw[islice]-0.5) < 1.e-12, f"Faces of the velocity along z not on the centreline: {zw[islice]}."
    data_ref = np.loadtxt("ghiav.txt") # https://gist.github.com/ivan-pi/caa6c6737d36a9140fbcf2ea59c78b3c
    data_interp = np.interp(data_ref[:,0],np.r_[0.,yp,1.],np.r_[0.,data[0,:,islice],0.])
    np.testing.assert_allclose(data_interp, data_ref[:,1], rtol=0., atol=err_max)
if __name__ == "__main__":
    test_ldc()
    print("Passed!")
