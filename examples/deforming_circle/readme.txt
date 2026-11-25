Deforming circle (2D advection)
===============================

A circle advected by a prescribed, time-dependent vortex that stretches it
into a thin filament and then reverses, so that the circle should return to
its initial shape after one period (R. J. LeVeque, SIAM J. Numer. Anal. 33
(1996) 627-665; W. J. Rider and D. B. Kothe, J. Comput. Phys. 141 (1998)
112-152).

Set-up
------
* Unit square in the y-z plane (one cell along x), 256 x 256 cells; circle
  of radius 0.15 centred at (0.5, 0.75).
* The flow is not solved (`is_solve_ns = F`): the velocity is that of
  `inivel = 'vor-x'`, with a period of 4.
* ACDI by default; other interface methods are selected with
  `ifc_reconstruction`.

Result
------
The phase indicator at t = 4 (`psi_fld_*` of the last step) against the
initial one: the L1 norm of the difference measures the error of the
advection. The test of the same name does this comparison.
