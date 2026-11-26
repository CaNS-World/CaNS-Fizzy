Rising bubble (2D)
==================

The two benchmark cases of S. Hysing et al., Int. J. Numer. Methods Fluids
60 (2009) 1259-1288: a bubble rising in a heavier fluid in a closed
container.

* `input-case-01.nml`: density ratio 10, viscosity ratio 10, surface tension
  24.5; the bubble ends up as an ellipsoidal cap.
* `input-case-02.nml`: density ratio 1000, viscosity ratio 100, surface
  tension 1.96; the bubble develops thin skirts.

Set-up
------
* Domain of 1 x 2 in the y-z plane (one cell along x), 64 x 128 cells; bubble
  of radius 0.25 centred at (0.5, 0.5); gravity 0.98 along -z; walls on all
  sides.
* ACDI by default; other interface methods are selected with
  `ifc_reconstruction`.
* Run with, e.g., `mpirun -n 4 ./cans input-case-01.nml`.

Result
------
At t = 3 the centre of mass of the bubble of case 1 is at about z = 1.08
(1.081 in the benchmark); its interface is compared with the reference
contour in the test `rising_bubble`. For case 2 the benchmark groups give
between 1.125 and 1.138.
