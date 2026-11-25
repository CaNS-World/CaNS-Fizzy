#!/usr/bin/env python
def test_contact_angle():
    import numpy as np
    #
    # static contact angle of 60 degrees at the lower wall, with the three interface reconstructions (THINC
    # with the AB2 and with the direction-split advection): a
    # two-dimensional drop that starts as the circular cap of that angle (the equilibrium shape without
    # gravity) must stay there over 2000 time steps (one time unit, about the capillary time scale of the
    # drop); with an angle of 90 degrees instead, the wetted width changes by 10 % and the velocity is
    # 0.05 in the same time
    #
    ny, nz = 128, 48
    dy, dz = 4./ny, 1.5/nz
    theta = np.radians(60.)
    area = np.pi*0.5**2/2. # (the cap has the area of half a disk of radius 0.5)
    r = np.sqrt(area/(theta-np.sin(theta)*np.cos(theta)))
    height_cap, width_cap = r*(1.-np.cos(theta)), 2.*r*np.sin(theta)
    for case in ['plic', 'thinc', 'thincsplit', 'acdi']:
        psi_beg = np.fromfile(f'data-{case}/psi_fld_0000000.bin').reshape((ny, nz), order='F')
        psi = np.fromfile(f'data-{case}/psi_fld_0002000.bin').reshape((ny, nz), order='F')
        assert np.all(np.isfinite(psi)), f"{case}: invalid values."
        height = 0.5*(np.sum(psi[ny//2-1,:])+np.sum(psi[ny//2,:]))*dz
        width = (1.5*np.sum(psi[:,0])-0.5*np.sum(psi[:,1]))*dy # extrapolated to the wall from the first two layers
        vel = max(np.max(np.abs(np.fromfile(f'data-{case}/ve{c}_fld_0002000.bin'))) for c in ['y', 'z'])
        assert abs(height/height_cap-1.) < 0.01, f"{case}: height of the drop {height}, expected {height_cap}."
        assert abs(width/width_cap-1.) < 0.01, f"{case}: wetted width {width}, expected {width_cap}."
        assert vel < 5.e-3, f"{case}: velocity too large: {vel} (below 1.e-3 is expected)."
        assert abs(np.sum(psi)/np.sum(psi_beg)-1.) < 1.e-12, f"{case}: volume not conserved."
if __name__ == "__main__":
    test_contact_angle()
    print("Passed!")
