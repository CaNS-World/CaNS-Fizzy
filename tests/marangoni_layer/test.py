#!/usr/bin/env python
def test_marangoni_layer():
    import numpy as np
    #
    # two superposed fluid layers between walls, with a flat interface at mid height and a surface tension that
    # varies with the scalar, sigma = sigma0 + ds*s, where s = cos(2*pi*y/l) does not diffuse and is hardly
    # advected: the tangential (Marangoni) force drives a steady cellular flow that, in the Stokes limit, has
    # the stream function g(z)*sin(2*pi*y/l), with g biharmonic in each layer; its eight coefficients follow
    # from no slip at the walls, no flow through the interface, a continuous tangential velocity, and the jump
    # of the shear stress, equal to the gradient of the surface tension
    #
    #   plic   : geometric reconstruction, equal viscosities, 16^2 cells
    #   plic-mu: the same with a viscosity ratio of 10, which the viscosity weighting of the force handles
    #   acdi   : diffuse interface, equal viscosities, 32^2 cells
    #
    # the tolerances are about 1.4 times the errors of the methods (second order for the geometric
    # reconstruction: 2.4 % with 16 cells, and first order for the diffuse interface: 9.8 % with 32 cells)
    #
    l = 1.; h = 1.
    cases = {"plic": (16,1.,1.,0.01,3.5e-2), "plic-mu": (16,1.,0.1,0.001,3.5e-2), "acdi": (32,1.,1.,0.01,1.4e-1)}
    def basis(z,w):
        ch = np.cosh(w*z); sh = np.sinh(w*z)
        return (np.array([ch,z*ch,sh,z*sh]),
                np.array([w*sh,ch+w*z*sh,w*ch,sh+w*z*ch]),
                np.array([w**2*ch,2*w*sh+w**2*z*ch,w**2*sh,2*w*ch+w**2*z*sh]))
    def exact(mu_top,mu_bot,ds,yv,zc,yc,zw):
        w = 2*np.pi/l
        a = np.zeros((8,8)); r = np.zeros(8) # coefficients of the lower (0:4) and upper (4:8) layers
        f,f1,f2 = basis(-h/2,w); a[0,0:4] = f; a[1,0:4] = f1
        f,f1,f2 = basis( h/2,w); a[2,4:8] = f; a[3,4:8] = f1
        f,f1,f2 = basis(0.,w)
        a[4,0:4] = f; a[5,4:8] = f
        a[6,0:4] = f1; a[6,4:8] = -f1
        a[7,4:8] = mu_top*f2; a[7,0:4] = -mu_bot*f2; r[7] = ds*w
        c = np.linalg.solve(a,r)
        def g(z,d):
            zz = z-h/2; out = np.zeros_like(zz)
            for layer,sl in ((zz < 0,slice(0,4)),(zz >= 0,slice(4,8))):
                out[layer] = np.tensordot(c[sl],basis(zz[layer],w)[d],axes=1)
            return out
        return np.outer(np.sin(w*yv),g(zc,1)),np.outer(-w*np.cos(w*yc),g(zw,0))
    for case,(n,mu1,mu2,ds,err_max) in cases.items():
        def load(field):
            data = np.fromfile("data-"+case+"/fld_"+field+".bin")
            assert np.all(np.isfinite(data)), f"{case}, {field}: invalid values."
            return data[:-2].reshape((1,n,n),order="F")[0] # the last two entries are the time and the time step number
        v = load("v"); w = load("w"); psi = load("psi")
        yc = (np.arange(n)+0.5)*l/n; yv = (np.arange(n)+1.)*l/n
        zc = (np.arange(n)+0.5)*h/n; zw = (np.arange(n)+1.)*h/n
        mu_top,mu_bot = (mu1,mu2) if psi[0,-1] > 0.5 else (mu2,mu1)
        v_ex,w_ex = exact(mu_top,mu_bot,ds,yv,zc,yc,zw)
        err = np.sqrt((np.mean((v-v_ex)**2)+np.mean((w-w_ex)**2))/(np.mean(v_ex**2)+np.mean(w_ex**2)))
        assert err < err_max, f"{case}: the error of the velocity is too large: {err} >= {err_max}."
if __name__ == "__main__":
    test_marangoni_layer()
    print("Passed!")
