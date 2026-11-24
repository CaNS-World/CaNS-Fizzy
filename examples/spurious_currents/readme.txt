Spurious currents of a circular drop (2D)
=========================================

A circular interface in equilibrium, for which the exact velocity is zero
(static cases) or uniform (moving case): the remaining velocity measures the
imbalance between the surface tension and the pressure gradient.

* `input.nml`: equal densities and viscosities, surface tension 3, fluids at
  rest, up to t = 5.
* `input-static.nml`: density ratio 1000, viscosity ratio 100, surface
  tension 0.05, fluids at rest, up to t = 0.2.
* `input-moving.nml`: as the previous one, with the drop carried across the
  periodic box by a uniform velocity along y (`inivel = 'uni-y'`).

Set-up
------
* Unit square in the y-z plane (one cell along x), periodic, 64 x 64 cells;
  circle of radius 0.2 at the centre.
* ACDI by default; other interface methods are selected with
  `ifc_reconstruction`.

Result
------
The largest velocity (minus the uniform one in the moving case), e.g., as a
capillary number mu*|u|/sigma, should stay small and not grow in time.
