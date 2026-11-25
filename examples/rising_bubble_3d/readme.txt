Rising bubble (3D)
==================

The three-dimensional counterparts of the two cases of the 2D rising bubble
(see `examples/rising_bubble_2d`), with the same fluid properties:
`input-case-01.nml` (density ratio 10) and `input-case-02.nml` (density
ratio 1000).

Set-up
------
* Box of 1 x 1 x 2 with walls on all sides, 128 x 128 x 256 cells; bubble of
  radius 0.25 centred at (0.5, 0.5, 0.5); gravity 0.98 along -z.
* ACDI by default; other interface methods are selected with
  `ifc_reconstruction`.
* Run with, e.g., `mpirun -n 4 ./cans input-case-01.nml`.

Result
------
The bubble rises and deforms until t = 3; its rise velocity and shape are
obtained from `psi_fld_*` and `vez_fld_*`.
