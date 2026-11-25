#!/usr/bin/env python
def test_rising_bubble():
    import numpy as np
    from read_single_field_binary import read_single_field_binary
    #
    # comparing rising bubble interface location at t=3 against
    # reference data from Hysing et al. IJNMF 60.11 (2009): 1259-1288;
    # the largest distance of the interface to the reference contour is between 0.005 and 0.009 for the
    # interface methods and pressure solvers of the code, on this grid of 64 x 128 cells
    #
    data,xp,yp,zp,xu,yv,zw = read_single_field_binary("data/","psi_fld_0006000.bin",np.array([1,1,1]))
    import matplotlib
    import matplotlib.pyplot as plt
    from scipy.spatial.distance import directed_hausdorff
    cs=plt.contour(yp,zp,np.transpose(data[0,:,:]),levels=[0.5])
    p = cs.get_paths()[0] if hasattr(cs,'get_paths') else cs.collections[0].get_paths()[0]
    v = p.vertices
    plt.clf()
    v_ref=np.loadtxt("c1g2l3s.txt")
    err = directed_hausdorff(v,v_ref)[0]
    err_max = 2.e-2
    assert err < err_max, f"Error is too large: {err} >= {err_max}."
    #
    # plotting
    #
    #x_psi = v[:,0]
    #y_psi = v[:,1]
    #x_psi_ref = v_ref[:,0]
    #y_psi_ref = v_ref[:,1]
    #plt.plot(x_psi    ,y_psi    )
    #plt.plot(x_psi_ref,y_psi_ref)
    #plt.show()
if __name__ == "__main__":
    test_rising_bubble()
