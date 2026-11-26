Lid-driven cavity (2D, single fluid)
====================================

Flow in a square cavity driven by its top wall, at a Reynolds number of 100
(U. Ghia, K. N. Ghia and C. T. Shin, J. Comput. Phys. 48 (1982) 387-411).

Set-up
------
* Unit square in the y-z plane (one cell along x), 64 x 64 cells; the wall
  at the top moves along y with unit velocity.
* A single fluid (`inipsi = 'uni'`) with viscosity 0.01.

Result
------
The flow is steady at the end of the run (t = 15.5). The velocity profiles
along the two centrelines agree with the tabulated ones of Ghia et al.; the
test of the same name does this comparison.
