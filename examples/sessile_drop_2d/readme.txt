Sessile drop with a contact angle (2D)
======================================

A drop on a wall relaxes from a half disc (a contact angle of 90 degrees) to
the circular cap of the prescribed static contact angle, 60 degrees, without
gravity.

Set-up
------
* Domain of 4 x 1.5 in the y-z plane (one cell along x), 256 x 96 cells;
  periodic along y, walls along z; half disc of radius 0.5 on the lower
  wall.
* `contact_angle(0:1,1:3) = 90.,90., 90.,90., 60.,90.`: 60 degrees on the
  lower wall.
* Geometric reconstruction (`ifc_reconstruction = 'plic'`); the other methods
  take the angle too.

Result
------
The drop spreads until it is a circular cap of the same area that meets the
wall at 60 degrees: height 0.400 and wetted width 1.385, for the initial
radius of 0.5. The test `contact_angle` starts from that cap and checks that
it stays in equilibrium.
