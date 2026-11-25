this is the directory where the tests are stored; they are run on GitHub for each pull request to `main` and
each push to it (see `.github/workflows/`)

each test is run from its own directory with `bash testit.sh`, and uses the executable `run/cans` built beforehand
(except `plic_unit`, which builds its own small program, and `single_precision`, which also builds the
code in single precision, in a copy of the source tree, and compares a few cases with double precision);
`testall.sh` builds the variants of the code that the tests need (the default build, the one without the
balanced capillary pressure split, and the one with the HYPRE solver) and runs the tests of each one, those that
are shared by several interface methods once per method (the method is an input: the variants `acdi`, `thinc`,
`thincsplit` and `plic` use the same build)

what the tests need: `mpirun` with 4 tasks, and Python with `numpy`, `scipy`, `matplotlib` and `pytest`;
they run in `run/`, so one at a time; they assume a double-precision build (the output is read as 8-byte reals,
and several tolerances are at double-precision round-off): only `poisson` and `single_precision` apply to a
single-precision one

environment variables:
  TEST_TWO_FLUID_INPUT: a line added to the `two_fluid` namelist of the tests shared by several interface
                        methods, e.g., "ifc_reconstruction = 'plic'" (set by `testall.sh` and by the CI)
  TEST_SHORT          : if set, the three long cases (`lid_driven_cavity`, `differentially_heated_cavity`,
                        `deforming_circle`) run a fraction of their time steps without their comparison, and
                        `single_precision` is skipped (the debug builds of the CI)
  see `testall.sh` for those of that script (BUILD_FLAGS, TESTS, NJOBS)

the tests, and what each checks:
  lid_driven_cavity            : centreline velocity profiles against Ghia et al. (single fluid)
  differentially_heated_cavity : Nusselt number at the wall against its reference value (single fluid, scalar)
  taylor_green_vortex_2d       : orders of convergence in space and time, in the three orientations
  poisson                      : residual of the pressure solver, for 79 combinations of grid, tasks, pencils and
                                 boundary conditions
  start_stop                   : a run against the same run stopped and restarted (single fluid)
  navier_slip                  : channel flow with a slip length against the exact profile
  forced_hit                   : the two large-scale forcings, from rest, against their solutions
  zero_net_acceleration        : balances of the forcing that keeps the mean momentum, and rejected inputs
  forcing_targets              : forcing that keeps the mean of a scalar, and the target bulk velocity and scalar
  deforming_circle             : a circle in a vortex that reverses returns to its initial shape
  rising_bubble                : interface of a rising bubble against the benchmark of Hysing et al., case 1
  restart_two_fluid            : a two-fluid run with three scalars against the same run restarted, bitwise
  implicit_diffusion           : implicit diffusion along z and in all directions against explicit runs
  markers                      : two markers reproduce one, and drops of different markers do not merge
  marangoni_layer              : thermocapillary flow of two layers against the Stokes solution
  marangoni_markers            : thermocapillary drops in several markers reproduce those in one
  plic_unit                    : kernels of the geometric reconstruction and of THINC in the sweeps
  plic_advection               : reversed advection in 2D and 3D: shape error, volume, bounds
  plic_momentum                : a drop with a density ratio of 1e6 in a uniform flow; restart
  viscosity_limiter            : the start of the Hysing case 2 stays bounded with the geometric method
  contact_angle                : a drop with a contact angle of 60 degrees stays in equilibrium
  options                      : options that no other test selects, on a short rising bubble
  single_fluid                 : a single fluid set as phase 1, with the interface not tracked, against phase 2
  input_checks                 : the default input file runs, and invalid inputs are rejected with a message
  single_precision             : eight short cases in single against double precision
