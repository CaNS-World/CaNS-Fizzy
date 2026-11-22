Thermocapillary migration of a drop (2D)
==========================================

A drop in a fluid at rest with a uniform temperature gradient along z. The
surface tension decreases with the temperature (`dsigma_ds < 0`), and the
drop migrates towards the hot wall.

Case of C. Ma and D. Bothe, Int. J. Multiphase Flow 37 (2011) 1045-1058, also
used by I. Seric, S. Afkhami and L. Kondic, J. Comput. Phys. (2018),
arXiv:1703.00327: Reynolds and Marangoni numbers of 0.72 and capillary number of
0.0576; the drop has half the density, viscosity, conductivity and heat
capacity of the ambient fluid.

Set-up
------
* Units: drop radius R, density and viscosity of the ambient fluid, and the
  velocity U = |dsigma_ds|*G*R/mu, with G the temperature gradient.
* Box of 4 x 4 radii, drop at its centre; walls with prescribed temperature
  at the bottom (0) and at the top (4), free slip and zero heat flux at the
  sides (periodic along x, with one cell).
* `s_ref = 2` is the temperature at the initial position of the drop, where
  the surface tension is `sigma`.
* 16 grid cells per radius, geometric reconstruction; other interface
  methods are selected with `ifc_reconstruction`.

Result
------
The drop velocity (the volume average of the vertical velocity over the
drop), in units of U, rises and levels off after about one time unit at
about 0.147 (the same with 16 and 32 cells per radius).
The theory of Young, Goldstein and Block (J. Fluid Mech. 6 (1959) 350-356)
for an unbounded fluid at zero Reynolds and Marangoni numbers gives 0.229;
the lower value is due to the small box and to the finite Reynolds and
Marangoni numbers.
