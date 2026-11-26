Drop or bubble in a pressure-driven channel flow (3D)
=====================================================

A single spherical drop (`input-droplet_laden.nml`, ten times denser than
the carrier fluid) or bubble (`input-bubble_laden.nml`, ten times lighter)
in a plane channel driven by a constant body force.

Set-up
------
* Channel of 2 x 2 x 1, periodic along x and y, walls along z, 128 x 128 x 64
  cells; drop or bubble of radius 0.25 at the centre.
* Flow at rest initially, driven by `bforce(1) = 1`; equal viscosities, and
  a surface tension chosen for a capillary number of 0.1 (see the comments
  in the input files).
* Run with, e.g., `mpirun -n 4 ./cans input-droplet_laden.nml`.

Result
------
A set-up to start from for particle-laden channel flows, more than a
benchmark: the drop or bubble is carried and deformed by the developing
flow. More bodies are added as lines of `spheres.in`.
