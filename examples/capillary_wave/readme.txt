Capillary wave (2D)
===================

Damped oscillation of a small-amplitude wave on the interface between two
viscous fluids at rest, without gravity.

Case with the analytical solution of A. Prosperetti, Phys. Fluids 24 (1981)
1217-1223, for two fluids of the same kinematic viscosity: here the density
and viscosity ratios are both 10 (kinematic viscosity 0.01).

Set-up
------
* Domain of 1 x 3 in the y-z plane (one cell along x), 128 x 384 cells;
  periodic along y, walls along z.
* Initial condition `cap-wav-1d-x`: one wavelength along y with an amplitude
  of 1 % of it, interface at mid height, fluids at rest; surface tension 1.
* ACDI by default, with `psi_thickness_factor = 1`; other interface methods
  are selected with `ifc_reconstruction`.

Result
------
The amplitude of the wave, e.g., from the position of the interface at one
end of the domain, oscillates and decays in time (up to t = 10), to be
compared with the analytical solution.
