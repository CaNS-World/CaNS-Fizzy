Translating circle (2D advection)
=================================

A circle carried by a uniform velocity across a periodic box, back to its
initial position.

Set-up
------
* Unit square in the y-z plane (one cell along x), periodic, 256 x 256 cells;
  circle of radius 0.25 at the centre.
* The flow is not solved (`is_solve_ns = F`): the velocity is uniform and
  equal to one along y (`inivel = 'uni-y'`), so the circle crosses the box
  once per time unit, five times in the run.
* ACDI by default; other interface methods are selected with
  `ifc_reconstruction`.

Result
------
The phase indicator at t = 5 against the initial one measures the error of
the advection.
