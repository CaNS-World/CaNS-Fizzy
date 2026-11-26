Differentially heated cavity (2D, single fluid)
===============================================

Natural convection of air in a square cavity with a hot and a cold side
wall and adiabatic top and bottom walls, under the Boussinesq approximation:
Rayleigh number 1e6, Prandtl number 0.71 (G. de Vahl Davis, Int. J. Numer.
Methods Fluids 3 (1983) 249-264).

Set-up
------
* Unit square in the y-z plane (one cell along x), 128 x 128 cells; gravity
  along -z; temperature -0.5 and 0.5 at the two walls along y.
* Units: cavity size, free-fall velocity and temperature difference, so that
  the viscosity is sqrt(Pr/Ra) and the conductivity 1/sqrt(Pr*Ra).
* A single fluid (`inipsi = 'uni'`, equal properties), one scalar with
  `is_boussinesq_buoyancy = T`.

Result
------
The flow reaches a steady state within the 50000 time steps. The mean
Nusselt number at the cold wall is about 8.8 (the test of the same name
compares it with 8.8252, within 1 %).
