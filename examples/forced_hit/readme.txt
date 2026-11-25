Forced homogeneous isotropic turbulence (3D, single fluid)
=========================================================

A triply-periodic box with a steady large-scale forcing of Taylor-Green
type (`is_forced_hit = T`, see `src/forcing.f90`), which drives a flow from
rest.

Set-up
------
* Box of side 2 pi, 256 x 256 x 256 cells, periodic in all directions.
* A single fluid (`inipsi = 'uni'`) with viscosity 0.002.
* A dispersed phase can be added with `inipsi` and `spheres.in`.

Result
------
The flow develops under the forcing; the velocity fields are in
`ve?_fld_*`.
