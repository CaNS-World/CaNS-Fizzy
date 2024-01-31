#!/usr/bin/env python
def test_ldc():
    import numpy as np
    from read_single_field_binary import read_single_field_binary
    #
    # testing u
    #
    data,xp,yp,zp,xu,yv,zw = read_single_field_binary("vey_fld_0020000.bin",np.array([1,1,1]))
    islice = np.size(data[0,0,:])//2
    data_ref = np.loadtxt("ghiau.txt") # https://gist.github.com/ivan-pi/3e9326d18a366ffe6a8e5bfda6353219
    data_interp = np.interp(data_ref[:,0],zp,0.5*(data[0,islice,:]+data[0,islice+1,:]))
    np.testing.assert_allclose(data_interp, data_ref[:,1], rtol=1.e-1, atol=1.e-1)
    #
    # plotting
    #
    #import matplotlib
    #import matplotlib.pyplot as plt
    #plt.plot(data_ref[:,0],data_ref[:,1],'.')
    #plt.plot(zp,data[0,islice,:],'--k')
    #plt.show()
    #
    # testing v
    #
    data,xp,yp,zp,xu,yv,zw = read_single_field_binary("vez_fld_0020000.bin",np.array([1,1,1]))
    islice = np.size(data[0,:,0])//2
    data_ref = np.loadtxt("ghiav.txt") # https://gist.github.com/ivan-pi/caa6c6737d36a9140fbcf2ea59c78b3c
    data_interp = np.interp(data_ref[:,0],yp,0.5*(data[0,:,islice]+data[0,:,islice+1]))
    np.testing.assert_allclose(data_interp, data_ref[:,1], rtol=1.e-1, atol=1.e-1)
    #import matplotlib
    #import matplotlib.pyplot as plt
    #plt.plot(data_ref[:,0],data_ref[:,1],'.')
    #plt.plot(yp,data[0,:,islice],'--k')
    #plt.show()
if __name__ == "__main__":
    test_ldc()
