Zalesak disk (2D advection)
===========================

Solid-body rotation of a slotted disk (S. T. Zalesak, J. Comput. Phys. 31
(1979) 335-362): the slot and its sharp corners test how well the interface
method keeps the shape.

Set-up
------
* Unit square in the y-z plane (one cell along x), 256 x 256 cells; slotted
  disk of radius 0.15 (`inipsi = 'zalesak-disk-x'`).
* The flow is not solved (`is_solve_ns = F`): the velocity is a rotation
  about the centre of the domain with a period of one (`inivel = 'rot-x'`).
* ACDI by default; other interface methods are selected with
  `ifc_reconstruction`.

Result
------
The phase indicator after a whole number of turns (the run ends at t = 5,
after five) against the initial one measures the error of the advection.
