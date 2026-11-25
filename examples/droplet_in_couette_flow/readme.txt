Drop in a plane Couette flow (3D)
=================================

Deformation of a neutrally buoyant drop in a simple shear flow between two
walls moving in opposite directions.

Set-up
------
* Channel of 2 x 2 x 1, periodic along x and y, walls along z with velocities
  0.5 and -0.5 along x (unit shear rate), 128 x 128 x 64 cells; drop of
  radius 0.25 at the centre.
* Equal densities and viscosities; Reynolds number 20 and capillary number
  0.1, based on the drop radius and the shear rate (see the comments in the
  input file).
* The flow starts from the Couette profile (`inivel = 'cou'`).

Result
------
The drop deforms into an ellipsoid inclined to the flow and reaches a steady
shape; its deformation parameter, (L-B)/(L+B) with L and B the longest and
shortest axes in the x-z plane, can be measured from `psi_fld_*`.
