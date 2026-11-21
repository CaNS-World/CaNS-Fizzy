# About the input file `input.nml`

Consider the following input file as an example (corresponding to a turbulent plane channel flow with bubbles, and walls at different temperatures). `&dns` defines a Fortran namelist containing the necessary computational parameters to set a case.


```fortran
&dns
ng(1:3) = 512, 256, 144
l(1:3) = 6., 3., 1.
gtype = 1, gr = 0.
cfl = 0.95, dtmax = 1.e5, dt_f = -1.
is_solve_ns = T, is_track_interface = T
inivel = 'poi'
is_wallturb = T, is_forced_hit = F
nstep = 100000, time_max = 100., tw_max = 0.1
stop_type(1:3) = T, F, F
restart = F, is_overwrite_save = T, nsaves_max = 0
icheck = 10, iout0d = 10, iout1d = 100, iout2d = 500, iout3d = 10000, isave = 5000
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,2) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,3) =  0.,0.,   0.,0.,   0.,0.
bcpre(0:1,1:3  ) =  0.,0.,   0.,0.,   0.,0.
bforce(1:3) = 0., 0., 0.
is_forced(1:3) = T, F, F
gacc(1:3) = 0., 0., 0.
nscal = 1
dims(1:2) = 2, 2, ipencil_axis = 1
/
```

In addition, the `&two_fluid` namelist contains the physical and transport properties of the two phases and the parameters of the interface-capturing method, and the `&scalar` namelist, read when `nscal = 1`, the parameters of the transport of a scalar (e.g., temperature):

```fortran
&scalar
inisca          = 'cou'
cbcsca(0:1,1:3) = 'P','P',  'P','P',  'D','D'
bcsca(0:1,1:3)  =  0.,0. ,   0.,0. ,   1.,0.
ssource         = 0.
is_boussinesq_buoyancy = F
/

&two_fluid
inipsi            = 'bub3'
cbcpsi(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
cbcnor(0:1,1:3,1) = 'P','P',  'P','P',  'N','N'
cbcnor(0:1,1:3,2) = 'P','P',  'P','P',  'N','N'
cbcnor(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
bcpsi(0:1,1:3)    =  0.,0. ,   0.,0. ,   0.,0.
bcnor(0:1,1:3,1)  =  0.,0.,    0.,0.,    0.,0.
bcnor(0:1,1:3,2)  =  0.,0.,    0.,0.,    0.,0.
bcnor(0:1,1:3,3)  =  0.,0.,    0.,0.,    0.,0.
sigma             =  1.
rho12(1:2)        =  1., 1.
mu12(1:2)         =  0.01, 0.01
ka12(1:2)         =  0.01, 0.01
cp12(1:2)         =  1., 1.
beta12(1:2)       =  0., 0.
psi_thickness_factor = 0.51
/
```
<details>

<summary>Tip for vim/nvim users</summary>
Consider adding the following lines in your `.vimrc` file for syntax highlighting of the namelist file:

```vim
if has("autocmd")
  au BufNewFile,BufRead *.nml set filetype=fortran
  au BufNewFile,BufRead *.namelist set filetype=fortran
endif
```

</details>

---
---

```fortran
ng(1:3) = 512, 256, 144
l(1:3) = 6., 3., 1.
gtype = 1, gr = 0.
```

These lines set the computational grid.

`ng(1:3)` and `l(1:3)` are the **number of points**  and **domain length** in each direction.

For 2D simulations, use one grid cell with periodic boundary conditions for all variables along the homogeneous direction. Initialize the velocity normal to that direction to zero and leave its forcing disabled. Two periodic cells are supported as well.

The admissible MPI process grids depend on the homogeneous direction:

| Homogeneous direction | Grid | `dims(1:2)` | Process-count limit |
|---|---|---|---|
| X | `1, Ny, Nz` | `1, P` | `P <= min(Ny,Nz)` |
| Y | `Nx, 1, Nz` | `1, 1` | one MPI rank |
| Z | `Nx, Ny, 1` | `P, 1` | `P <= min(Nx,Ny)` |

`dims(1:2) = 0, 0` selects these process grids automatically.

With `ng(3)=1`, use `is_poisson_dtdma = F`, and use a uniform Z grid.

Grid and timestep formulas remain three-dimensional. Since the spacing in a one-cell direction equals its domain length, choose that length at least as large as the largest spacing in the two retained directions. Together with zero normal velocity, this avoids an artificial timestep restriction from the homogeneous direction. See the supplied 2D tests and examples, which use one cell along X; their two-dimensional initial conditions take the suffix `-x` (see `inivel`, `inisca` and `inipsi` below).

`gtype` and `gr` are the **grid stretching type** and **grid stretching parameter** that tweak the non-uniform grid in the third direction; zero `gr` implies no stretching. See `initgrid.f90` for more details. The following options are available for `gtype`:

* `1`: grid clustered towards both ends (default)
* `2`: grid clustered towards the lower end
* `3`: grid clustered towards the upper end
* `4`: grid clustered towards the middle

---

```fortran
cfl = 0.95, dtmax = 1.e5, dt_f = -1.
```

This line controls the simulation time step size.

The time step size is set to be equal to `min(cfl*dt_cfl,dtmax)` if `dt_f < 0`, and to `dt_f` otherwise. In the former case, the code prescribes the minimum value between `dtmax` and `cfl` times the maximum allowable time step `dt_cfl` (computed every `icheck` time steps; see below). `dtmax` is therefore used when a constant time step, smaller than `cfl*dt_cfl`, is required. If not, it should be set to a high value so that the time step is dynamically adjusted to `cfl*dt_cfl`. Finally, a constant time step size may be forced, irrespective of the temporal stability evaluation through `dt_f`.

Note: the diffusion terms limit `dt_cfl` to `1/(4*nu*(1/dx**2+1/dy**2+1/dz**2))`, the stability limit of the AB2 scheme, with the largest diffusivity `nu` of the two fluids and the smallest grid spacing along each direction. The two cells of the invariant direction of a two-dimensional run count too, so a large domain length along that direction avoids restricting the time step.

---

```fortran
is_solve_ns = T, is_track_interface = T
```

`is_solve_ns`, if true, enables the solution of the Navier-Stokes momentum equation. It can be set to false when a prescribed velocity field is desired.

`is_track_interface`, if true, enables the advection of the two-fluid interface.

The viscosities and densities of the two fluids are set in the `&two_fluid` namelist (see below).

---

```fortran
inivel = 'poi'
is_wallturb = T, is_forced_hit = F
```

These lines set the initial velocity field.

`inivel` **chooses the initial velocity field**. The following options are available:

* `zer`: zero velocity field
* `uni`: uniform velocity field equal to `uref`
* `cou`: plane Couette flow profile using the prescribed velocities at the normal walls
* `poi`: plane/duct Poiseuille flow profile with mean velocity `uref`
* `tbl`: temporal boundary layer profile with wall velocity `uref`
* `pdc`: plane Poiseuille flow profile with constant pressure gradient
* `log`: logarithmic channel/duct profile with mean velocity `uref`
* `hcp`: half channel with plane Poiseuille profile and mean velocity `uref`
* `hcl`: half channel with logarithmic profile and mean velocity `uref`
* `hdc`: half plane Poiseuille flow profile with constant pressure gradient
* `tgv`: three-dimensional Taylor-Green vortex
* `tgv-2d-x`: two-dimensional Taylor-Green vortex in the Y–Z plane, with zero X velocity
* `tgv-2d-y`: two-dimensional Taylor-Green vortex in the Z–X plane, with zero Y velocity
* `tgv-2d-z`: two-dimensional Taylor-Green vortex in the X–Y plane, with zero Z velocity
* `ant`: three-dimensional Antuono vortex
* `rot-x`, `rot-y`, `rot-z`: solid-body rotation with angular velocity $2 \pi$ around the axis X, Y or Z
* `vor-x`, `vor-y`, `vor-z`: two-dimensional, time-dependent serpentine vortex in the plane normal to X, Y or Z
* `def`: three-dimensional, time-dependent deformation field of LeVeque, with a period of 3

The profile initializers (`zer`, `uni`, `cou`, `poi`, `tbl`, `iop`, `log`, `hcl`, `hcp`, `pdc`, `hdc`) accept `-x`, `-y`, or `-z` suffixes, for example `poi-y`. A bare name always means `-x`.

| Suffix | Flow direction | Main profile variation | Other transverse direction |
|---|---|---|---|
| `-x` | X | Z | Y |
| `-y` | Y | Z | X |
| `-z` | Z | X | Y |

Wall velocities and forcing are read from the corresponding physical entries of `bcvel`, `bforce` and `is_forced`; these inputs and the grid are not reordered. The reference velocity `uref` is equal to one.

`is_wallturb`, if true, **superimposes a high amplitude disturbance on the initial velocity field** that effectively triggers transition to turbulence in a wall-bounded shear flow.

See `initflow.f90` for more details. The subroutines in this file may be easily extended to accommodate other initial conditions.

`is_forced_hit`, if true, introduces a homogeneous and isotropic forcing in the domain. See `forcing.f90` for more details.

---

```fortran
nstep = 100000, time_max = 100., tw_max = 0.1
stop_type(1:3) = T, F, F
restart = F, is_overwrite_save = T, nsaves_max = 0
```

These lines set the simulation termination criteria and whether the simulation should be restarted from a checkpoint file.

`nstep` is the **total number of time steps**.

`time_max` is the **maximum physical time**.

`tw_max` is the **maximum total simulation wall-clock time**.

`stop_type` sets which criteria for terminating the simulation are to be used (more than one can be selected, and at least one of them must be `T`)

* `stop_type(1)`, if true (`T`), the simulation will terminate after `nstep` time steps have been simulated
* `stop_type(2)`, if true (`T`), the simulation will terminate after `time_max` physical time units have been reached
* `stop_type(3)`, if true (`T`), the simulation will terminate after `tw_max` simulation wall-clock time (in hours) has been reached

`restart`, if true, **restarts the simulation** from a previously saved checkpoint, made of one file per field, `fld_<var>.bin` (e.g., `fld_u.bin`, `fld_psi.bin`; `.h5` / `.bp` for the other supported backends), all from the same save.

`is_overwrite_save`, if true, overwrites the checkpoint files `fld_<var>.bin` at every save; if false, each checkpoint is saved to files `fld_???????_<var>.bin` (with `???????` denoting the time step number), and symbolic links `fld_<var>.bin` point to the last ones; to restart from another checkpoint, point the links to its files (e.g., `ln -sf fld_0000100_u.bin fld_u.bin`).

`nsaves_max`, if positive and `is_overwrite_save` is true, keeps instead the last `nsaves_max` checkpoints in rotating slots, `fld_????_<var>.bin` (with `????` denoting the slot number), with links `fld_<var>.bin` to the last one, and a log of the saves in `log_checkpoints.out`; a restarted run continues the rotation.

Finally, note that a checkpoint file will be saved before the simulation ends. Checkpoints are written to temporary files (e.g., `fld_u.tmp.bin`), renamed once complete, so that an interrupted save leaves the previous checkpoint intact (except with the ADIOS2 backend). See `is_full_restart` under `&io` for restarting the time integration exactly.

---

```fortran
icheck = 10, iout0d = 10, iout1d = 100, iout2d = 500, iout3d = 10000, isave = 5000
```

These lines set the frequency of time step checking and output:

* every `icheck` time steps **the new time step size** is computed according to the stability criterion and `cfl` (above), and the **velocity divergence is checked**: if it is too large or NaN, the simulation is aborted without saving a checkpoint (see also `is_mask_divergence_check` under `&other_options`)
* every `iout0d` time steps **history files with global scalar variables** are appended; e.g., the forcing and its balance, and the time step history
* every `iout1d` time steps **1d profiles** are written (e.g. velocity and its moments) to a file
* every `iout2d` time steps **2d slices of a 3d scalar field** are written to a file
* every `iout3d` time steps **3d scalar fields** are written to a file
* every `isave`  time steps a **checkpoint file** is written (`fld_???????.bin`), and a symbolic link for the restart file, `fld.bin`, will point to this last save so that, by default, the last saved checkpoint file is used to restart the simulation

1d, 2d, and 3d outputs can be tweaked by modifying files `out?d.h90` and recompiling the source. See also `output.f90` for more details. _Set any of these variables to `0` to skip the corresponding operation._

```fortran
tout0d = 0., tout1d = 0., tout2d = 0., tout3d = 0., tsave = 0.
```

These optional entries of `&dns` set the same outputs in **physical time**: a positive period (default `0`) takes priority over the corresponding time-step interval above. An output is written at the end of the first time step that reaches a multiple of the period, at the actual time of that step (the time step is not adjusted).

---

```fortran
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,2) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,3) =  0.,0.,   0.,0.,   0.,0.
bcpre(0:1,1:3  ) =  0.,0.,   0.,0.,   0.,0.
```

These lines set the boundary conditions (BC).

The **type** (BC) for each field variable is set by a row of six characters, `X0 X1  Y0 Y1  Z0 Z1`, where:

* `X0` `X1` set the type of BC for the field variable at the **lower** and **upper** boundaries in `x`
* `Y0` `Y1` set the type of BC for the field variable at the **lower** and **upper** boundaries in `y`
* `Z0` `Z1` set the type of BC for the field variable at the **lower** and **upper** boundaries in `z`

The four rows correspond to the three velocity components, and pressure, i.e. `u`, `v`, `w`, and `p`. The boundary conditions of the scalar and of the interface fields are set in the `&scalar` and `&two_fluid` namelists.

The following options are available:

* `P` periodic
* `D` Dirichlet
* `N` Neumann

The **last four rows** follow the same logic, but now for the BC **values** (dummy for a periodic direction).

---

```fortran
bforce(1:3) = 0., 0., 0.
is_forced(1:3) = T, F, F
gacc(1:3) = 0., 0., 0.
```

These lines set the flow forcing.

`bforce` is a constant **body force density term** in the direction in question (e.g., the negative of a constant pressure gradient) that can be added to the right-hand side of the momentum equation. The three values correspond to the three domain directions.

`is_forced`, when true in a direction, adds a dynamically computed uniform correction to `bforce`. The controller balances the prescribed force, the AB2-extrapolated viscous wall traction, and the discrete integrals of gravity, surface tension, and Boussinesq buoyancy. It preserves the initial domain-averaged mixture momentum, equivalently the streamwise-averaged mass-flux density, when the remaining pressure and conservative terms have zero global integral. It does not prescribe a target velocity or necessarily preserve volume flux or the instantaneous mass flux through every cross-section.

A dynamically forced direction must be periodic. Its transverse boundaries must be periodic or impermeable, since momentum flux through an open boundary is not included in the balance.

`gacc` is a constant **acceleration** in the direction in question (e.g., gravity) added to the right-hand-side of the momentum equation. The three values correspond to three domain directions.

When fixed or dynamic forcing is active, `forcing.out` contains time, `bforce`, the three equivalent pressure gradients from the dynamic correction, and the three volume-averaged velocities. With dynamic forcing, `forcing_balance.out` additionally contains time followed by six three-component vectors: wall compensation, gravity, surface tension, Boussinesq buoyancy, total applied force density, and the measured change in total mixture momentum per unit time. The last vector is therefore a final residual that also exposes contributions from terms assumed to have zero global integral, such as pressure and conservative advection.

At unequal density, the constant-coefficient pressure split does not in general have a zero global mixture-momentum contribution. In that configuration the controller still closes the listed source-force balance, while the last vector in `forcing_balance.out` reports the remaining pressure contribution. The variable-coefficient pressure solve does not have this limitation.

---

```fortran
dims(1:2) = 2, 2, ipencil_axis = 1
```

This line sets the domain decomposition and orientation of the computational subdomains.

`dims` is the **processor grid**, the number of domain partitions along the first and second decomposed directions (which depend on the selected default pencil orientation). `dims(1)*dims(2)` corresponds therefore to the total number of computational subdomains. Setting `dims(:) = [0,0]` will trigger a runtime autotuning step to find the processor grid that minimizes transpose times. Note, however, that other components of the algorithm (e.g., collective I/O) may also be affected by the choice of processor grid.

`ipencil_axis` sets the **orientation of the computational subdomains** (or pencils), being one of `[1,2,3]` for `[X,Y,Z]`-aligned pencils. X-aligned is the default if this option is not set, and should be optimal for all cases except for Z-implicit diffusion, where Z-pencils are recommended if `dims(2) > 1` in the input file; see the description of the `&numerics` namelist below.

# About the `&cudecomp` namelist under `input.nml`

In addition to the `&dns` namelist in the input file, there is an **optional** namelist to set some runtime configurations for the *cuDecomp* library. Consider the following `&cudecomp` namelist, which corresponds to the default options in case the file is not provided:

```fortran
&cudecomp
cudecomp_t_comm_backend = 0, cudecomp_is_t_enable_nccl = T, cudecomp_is_t_enable_nvshmem = T
cudecomp_h_comm_backend = 0, cudecomp_is_h_enable_nccl = T, cudecomp_is_h_enable_nvshmem = T
```

The first line sets the configuration for the transpose communication backend autotuning. Here `cudecomp_t_comm_backend` can be one of:

* `1` -> `CUDECOMP_TRANSPOSE_COMM_MPI_P2P`
* `2` -> `CUDECOMP_TRANSPOSE_COMM_MPI_P2P_PL`
* `3` -> `CUDECOMP_TRANSPOSE_COMM_MPI_A2A`
* `4` -> `CUDECOMP_TRANSPOSE_COMM_NCCL`
* `5` -> `CUDECOMP_TRANSPOSE_COMM_NCCL_PL`
* `6` -> `CUDECOMP_TRANSPOSE_COMM_NVSHMEM`
* `7` -> `CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL`
* `8` -> `CUDECOMP_TRANSPOSE_COMM_NVSHMEM_SM`
* any other value -> enable runtime transpose backend autotuning

The other two boolean values, enable/disable the NCCL (`cudecomp_is_t_enable_nccl`) and NVSHMEM (`cudecomp_is_t_enable_nvshmem`) options for *transpose* communication backend autotuning.

The second line is analogous to the first one, but for halo communication backend autotuning. Here `cudecomp_h_comm_backend` can be one of:

* `1` -> `CUDECOMP_HALO_COMM_MPI`
* `2` -> `CUDECOMP_HALO_COMM_MPI_BLOCKING`
* `3` -> `CUDECOMP_HALO_COMM_NCCL`
* `4` -> `CUDECOMP_HALO_COMM_NVSHMEM`
* `5` -> `CUDECOMP_HALO_COMM_NVSHMEM_BLOCKING`
* any other value -> enable runtime halo backend autotuning

The other two boolean values, enable/disable the NCCL (`cudecomp_is_h_enable_nccl`) and NVSHMEM (`cudecomp_is_h_enable_nvshmem`) options for *halo* communication backend autotuning.

Finally, it is worth recalling that passing `dims(1:2) = [0,0]` under `&dns` will trigger the *processor grid* autotuning, so there is no need to provide that option in the `&cudecomp` namelist.

# About the `&scalar` namelist under `input.nml`

This namelist is **optional** and defines the parameters needed to solve the transport equation of a scalar (e.g., temperature). It is read when the scalar is enabled in the `&dns` namelist:

```fortran
&dns
! (...)
nscal = 1
/
```

The default value of `nscal` is zero, and at most one scalar is supported.

---

The following example namelist defines the parameters needed for the differentially heated cavity case that may be found under `examples/`.

```fortran
&scalar
inisca          = 'dhc-x'
cbcsca(0:1,1:3) = 'P','P',  'D','D',  'N','N'
bcsca(0:1,1:3)  =  0.,0.,   -0.5,0.5,   0.,0.
ssource         = 0.
is_boussinesq_buoyancy = T
/
```

---

```fortran
inisca          = 'dhc-x'
```

This line sets the **initial scalar field**. The following options are available:

* `zer`: uniform scalar field equal to zero
* `uni`: uniform scalar field equal to `sref`, the mean of the boundary values along Z
* `cou`: linearly varying scalar field from bottom (`z = 0`) to top   (`z = l(3)`)
* `dhc`: linearly varying scalar field from left   (`x = 0`) to right (`x = l(1)`); with the suffix `-x`, i.e., `dhc-x`, the field varies along Y instead, for a cavity in the Y-Z plane

See `initflow.f90` for more details.

---

```fortran
cbcsca(0:1,1:3) = 'P','P',  'D','D',  'N','N'
bcsca(0:1,1:3)  =  0.,0.,   -0.5,0.5,   0.,0.
```

These lines set the **boundary conditions** for the scalar, just like for the pressure.

---

```fortran
ssource         = 0.
```

`ssource` is a **uniform volumetric source term** added to the right-hand side of the scalar equation.

---

**Buoyancy effects**

```fortran
is_boussinesq_buoyancy = T
```

`is_boussinesq_buoyancy`, if true, **activates the buoyancy term** of the so-called _Boussinesq approximation_ in the momentum equation, with the scalar as the temperature difference. The thermal expansion coefficient of each fluid is set by `beta12` in the `&two_fluid` namelist, and the gravitational acceleration by `gacc` in the `&dns` namelist, e.g.:

```fortran
&dns
! (...)
gacc(1:3) = 0., 0., -1.
nscal = 1
/
```

which, in this example, defines a negative gravitational acceleration along `z` with unit magnitude. Its default value is `[0., 0., 0.]`.

# About the `&two_fluid` namelist under `input.nml`

This namelist defines the initial condition and the boundary conditions of the interface, the physical and transport properties of the two fluids, and the interface thickness.

```fortran
inipsi            = 'bub3'
```

`inipsi` **chooses the initial volume fraction field for phase 1** of the two-fluid system. The following options are available:

* `zer` : uniform (zero) volume fraction field
* `uni` : uniform (unitary) volume fraction field
* `bub1`: planar film (lighter phase)
* `bub2`: two-dimensional (lighter phase) bubbles
* `bub3`: three-dimensional (lighter phase) bubbles
* `drp1`: planar film (heavier phase)
* `drp2`: two-dimensional (heavier phase) droplets
* `drp3`: three-dimensional (heavier phase) droplets
* `dis1`: planar film (same as `drp1`)
* `dis2`: two-dimensional dispersed circles (same as `drp2`)
* `dis3`: three-dimensional dispersed spheres (same as `drp3`)
* `cap-wav-1d`: planar small-amplitude capillary wave
* `zalesak-disk`: two-dimensional (lighter phase) slotted disk

The two-dimensional fields (`bub2`, `drp2`, `dis2`, `cap-wav-1d` and `zalesak-disk`) accept the suffixes `-x`, `-y` or `-z`, the direction normal to their plane; e.g., `bub2-x` are bubbles in the Y-Z plane. A bare name means `-y`, the X-Z plane.

With the geometric volume-of-fluid method (`INTERFACE_CAPTURING_SPLIT_VOF`), these fields are sharp: each cell takes the exact volume fraction of phase 1, instead of a smooth profile.

See `two_fluid.f90` for more details. For the `bub[1-3]`, `drp[1-3]`, and `dis[1-3]` initial fields, the position and size of the films/bubbles/droplets can be specified by a `spheres.in` file, where each line corresponds to an individual spherical/circular/planar bubble/droplet, specifying the cartesian coordinates of its center of mass, and its radius. **Note:** _The recommended way to initialize a dispersed multiphase flow is to use `dis?`, with phase 1 being the dispersed phase, and phase 2 the continuous one. This is important because `inivel` above considers phase 2 as the continuous one._
For example, the following file introduces in the computational domain two bubbles/droplets, the first centered at `[x,y,z] = [0.3,1.,0.5]` with a radius of `0.25`, and the second centered at `[1.5,0.6,1.]` with radius `0.5`.

```fortran
0.3 1.  0.5 0.25
1.5 0.6 1.  0.5
```

In the case of a planar film (the one-dimensional limit of the sphere), the radius will represent the film thickness, and the film will be aligned along the X-Y plane.

Finally, if the `spheres.in` file is not specified, the solver will initialize a single sphere, located at the center of the domain, with a size of `0.25` times the smallest domain length, excluding homogeneous directions for 2D and 1D initializations.

---

```fortran
cbcpsi(0:1,1:3)   = 'P','P',  'P','P',  'P','P'
cbcnor(0:1,1:3,1) = 'P','P',  'P','P',  'P','P'
cbcnor(0:1,1:3,2) = 'P','P',  'P','P',  'P','P'
cbcnor(0:1,1:3,3) = 'P','P',  'P','P',  'P','P'
bcpsi(0:1,1:3)    =  0.,0. ,   0.,0. ,   0.,0.
bcnor(0:1,1:3,1)  =  0.,0.,    0.,0.,    0.,0.
bcnor(0:1,1:3,2)  =  0.,0.,    0.,0.,    0.,0.
bcnor(0:1,1:3,3)  =  0.,0.,    0.,0.,    0.,0.
```

These lines set the boundary conditions (BC) for the volume fraction of phase 1 (`cbcpsi`) and for the three cartesian components of the interface normal (`cbcnor`).

The **type** (BC) is set by a row of six characters, `X0 X1  Y0 Y1  Z0 Z1` where,

* `X0` `X1` set the type of BC on the **lower** and **upper** boundaries in `x`
* `Y0` `Y1` set the type of BC on the **lower** and **upper** boundaries in `y`
* `Z0` `Z1` set the type of BC on the **lower** and **upper** boundaries in `z`

The following options are available:

* `P` periodic
* `D` Dirichlet
* `N` Neumann

The corresponding BC **values** (dummy for a periodic direction) are set in `bcpsi` for the volume fraction of phase 1 and in `bcnor` for the components of the interface normal.
For consistency, the components of the interface normal should be periodic when the volume fraction boundary conditions are periodic.

Note: to conserve mass in the domain, a solid wall boundary requires a homogeneous Neumann BC for the volume fraction and a homogeneous Dirichlet BC for the wall-normal component of the interface normal. In this case, the BC for the wall-tangential components of the interface normal should be set to homogeneous Neumann.

---

```fortran
sigma           = 1.
rho12(1:2)      = 10., 1.
mu12(1:2)       = 1.e-2, 1.e-3
ka12(1:2)       = 0., 0.
cp12(1:2)       = 0., 0.
beta12(1:2)     = 0., 0.
```

These lines specify the physical and transport properties of the two-phase system.

* `sigma` sets the value of the surface tension coefficient at the two-phase interface.
* `rho12` sets the value of the density of each phase. (Note: the initial conditions assume `rho12(1)` as the densest phase)
* `mu12` sets the value of the dynamic viscosity of each phase.
* `ka12` sets the value of the thermal conductivity of each phase.
* `cp12` sets the value of the constant pressure specfic heat capacity of each phase.
* `beta12` sets the value of the thermal expansion coeffcient of each phase.

By default, viscosity and thermal conductivity are interpolated arithmetically between the two phases. Setting `VISC_HARMONIC_INTERPOLATION=1` at compilation time switches both properties to harmonic phase interpolation.

---

```fortran
psi_thickness_factor = 0.51
```

This value sets the interface thickness (the higher, the thicker). It corresponds to the non-dimensional epsilon parameter (`eps`) from the ACDI method which, for ACDI, should always be larger than `0.5`; for the THINC/QQ VoF method, a corresponding `beta` parameter is related related to this factor by `beta = 1./(2*eps)`. This parameter controls the interface thickness as given by the hyperbolic tangent function: `h = 0.5*(1.+tanh(r/(2.*eps))) = 0.5*(1.+tanh(r*beta))`, with `r` being the signed distance to the interface in computational space. The recommended default values are `psi_thickness_factor = 0.51` for ACDI (`eps = 0.51`), and `psi_thickness_factor = 0.5` to `0.25` (`beta = 1` to `2`), for the THINC/QQ VoF.

Note: see the compile-time options `INTERFACE_CAPTURING_VOF` and `INTERFACE_CAPTURING_SPLIT_VOF` in [`docs/INFO_COMPILING.md`](INFO_COMPILING.md) for how to switch to VoF for interface capturing. The geometric method of `INTERFACE_CAPTURING_SPLIT_VOF` keeps the interface sharp, and does not use this parameter; every `iout0d` time steps, it appends to `log_vof.out` the time step number, the time, the volume of phase 1 and its relative change since the start of the run, the accumulated volume added by clipping round-off errors, the extrema of the volume fraction before clipping, the CFL numbers along the three directions, the number of sub-steps of the last time step, and the numbers of interfacial cells whose curvature comes from height functions, from the average of neighbouring cells, and without an estimate.

# About the `&numerics` namelist under `input.nml`

This namelist defines parameters related to the numerical discretization and computational method. The values below are the default ones, in case the namelist is not specified in the input file.

```fortran
&numerics
is_impdiff = F, is_impdiff_1d = F
is_poisson_dtdma = F
/
```

`is_impdiff_1d`, if true, enables the **implicit temporal integration of the diffusion along Z**, with the Crank-Nicolson scheme: the time step is then no longer restricted by the grid spacing along Z in the viscous (and scalar diffusion) terms, which pays off in wall-bounded flows with a fine grid near the walls at low Reynolds number or high viscosity. The terms integrated implicitly are the diffusion of each velocity component along Z due to its own gradient, with the local viscosity and density of the two-fluid mixture, and that of the scalar; the other viscous terms and the regularization of the interface remain explicit. Each line along Z is solved with a tridiagonal algorithm, which requires that **the domain is not decomposed along Z**: use `ipencil_axis = 3`, or `dims(2) = 1`; periodic boundaries along Z are not supported yet. The pressure field is then the one of the projection, and differs from the physical pressure by a small term proportional to the time step and the viscosity.

`is_impdiff`, if true and with `is_impdiff_1d = F`, enables the **implicit temporal integration of the diffusion in all directions**, so that the time step is not restricted by the diffusion terms. Its form follows the pressure solver of the build (see [`INFO_COMPILING.md`](INFO_COMPILING.md)):

* with the FFT-based solver (`CONSTANT_COEFFS_POISSON=1`, the default), **with constant coefficients**: the term integrated implicitly is the Laplacian of each velocity component (and of the scalar) times a reference diffusivity, the largest kinematic viscosity (thermal diffusivity) of the two fluids, and the remainder of the diffusion terms is integrated explicitly. The three (four, with a scalar) Helmholtz equations are solved with the FFT-based solver, as in *CaNS*, on CPUs and GPUs, with the same restrictions on the boundary conditions. This splitting is exact for a single-phase flow. With two fluids, the reference value must be the largest one for stability, so that the fluid with the smaller kinematic viscosity is integrated somewhat less accurately at time steps well beyond the explicit limit: in tests with rising bubbles at 10 times that limit, the errors were up to 4 times those of the integration with the local properties for a contrast of 100 between the kinematic viscosities, and similar for a contrast of 10.
* with the HYPRE solver (`CONSTANT_COEFFS_POISSON=0`), **with variable coefficients**: the terms integrated implicitly are the diffusion of each velocity component due to its own gradient, and that of the scalar, with the local properties of the mixture, as for `is_impdiff_1d = T` but in all directions; the other viscous terms remain explicit. The three (four) linear systems are solved with HYPRE, with the second entry of the inputs of the `&hypre` namelist (see below); their unknown is the change of the field due to the implicit terms, so that the tolerance is relative to the size of the diffusive update.

In both cases the interface regularization remains explicit.

Note: at time steps beyond the explicit limit of the diffusion terms, other errors may take over. In particular, the error of the constant-coefficient pressure splitting (`CONSTANT_COEFFS_POISSON=1`) grows with the time step and the density ratio: at a density ratio of 1000 it dominated the tests above, for all the implicit options, and was absent with the variable-coefficient pressure solver.

Note: `is_impdiff` is set to true automatically with `is_impdiff_1d = T`.

Finally, `is_poisson_dtdma`, if `.true.`, solves the Poisson/Helmholtz equations along Z with a distributed tridiagonal matrix algorithm (DTDMA). This approach may result in major gains in scalability for pencil-distributed simulations at scale, on many GPUs.

# About the `&hypre` namelist under `input.nml`

This namelist is **optional**, and tunes the iterative solver of the [HYPRE](https://github.com/hypre-space/hypre) library that is used for the variable-coefficient pressure Poisson equation, in builds with `CONSTANT_COEFFS_POISSON=0` (see [`INFO_COMPILING.md`](INFO_COMPILING.md)). Each input has two entries: the first for the pressure, and the second for the implicit diffusion in all directions (`is_impdiff = T`, see above). The values below are the default ones.

```fortran
&hypre
hypre_solver_i(1:2)  = 5, 5
hypre_tol(1:2)       = 1.e-6, 1.e-6
hypre_maxiter(1:2)   = 50, 50
hypre_precond_i(1:2) = 1, 1
hypre_precond_maxiter(1:2) = 1, 1
/
```

`hypre_solver_i` selects the **solver**:

* `1`: SMG, a geometric multigrid solver
* `2`: PFMG, a geometric multigrid solver with semicoarsening
* `3`: GMRES
* `4`: BiCGSTAB
* `5`: Flexible GMRES (default)

`hypre_tol` is the **relative tolerance** on the residual, and `hypre_maxiter` the **maximum number of iterations**; a warning is printed whenever the solver does not converge within it.

`hypre_precond_i` selects the **preconditioner** of the Krylov solvers (`3`, `4` and `5`), which is applied `hypre_precond_maxiter` times per iteration:

* `0`: SMG
* `1`: PFMG (default)
* `9`: none

The defaults were chosen from tests with three- and two-dimensional rising bubbles at density ratios up to 1000, with walls, periodic directions, stretched grids and an outflow. For a problem without a Dirichlet pressure boundary, the pressure is defined up to a constant, and the mean of the solution is set to zero; for such a problem, PFMG stops coarsening one level short of a single grid cell, where the operator would vanish and spoil the convergence. The same defaults suit the systems of the implicit diffusion, which are not singular and converge in a few iterations.

With an iterative solver, the **velocity divergence** follows the tolerance rather than the machine precision: with the default tolerance its maximum is typically `1.e-7` or smaller, and each factor of ten in `hypre_tol` lowers it by about as much, at a moderate cost (the solver gains roughly one digit per iteration). For this reason, the threshold test of the divergence check is skipped by default in these builds (`is_mask_divergence_check = T` under `&other_options`); it can be enforced again by setting it to false together with a tighter tolerance.

Note: a periodic direction must have an even number of grid cells (or a single one, for a two-dimensional case), a requirement of the multigrid solvers of HYPRE.

Note: other parameters (`hypre_gmres_k_dim`, `hypre_pfmg_relax_type`, `hypre_pfmg_num_pre_relax`, `hypre_pfmg_num_post_relax`, `hypre_pfmg_max_levels`) are meant for advanced use, and can be found in `param.f90`.

# About the `&other_options` namelist under `input.nml`

This namelist defines other parameters related to the monitoring, debugging, or benchmarking of a computation. The values below are the default ones, in case the namelist is not specified in the input file.

```fortran
&other_options
is_debug = T, is_timing = T
is_mask_divergence_check = F
/
```

In these lines, `is_debug` performs some **sanity checks for debugging purposes** that do not introduce computational overhead, and `is_timing` reports the wall-clock time per step.

Finally, `is_mask_divergence_check`, if true, skips the threshold test of the **velocity divergence check** performed every `icheck` time steps: by default, the simulation is aborted when the maximum divergence exceeds a small, precision-dependent threshold, which may be smaller than the round-off divergence of very large grids in single precision. A NaN in the velocity field always aborts the simulation, regardless of this option. In builds with the iterative pressure solver (`CONSTANT_COEFFS_POISSON=0`) the default is true (see the `&hypre` namelist).

Note: other parameters for `&numerics` and `&other_options` are not exposed here, as they are meant for very specific, advanced use or developers. They can be found in `param.f90`.

This optional namelist selects the **I/O backend** used for checkpointing and field output.

# About the `&io` namelist under `input.nml`

```fortran
&io
io_backend = 'mpiio'
/
```

The following options are currently available for `io_backend`:

* `mpiio`: raw binary output based on MPI-I/O (default)
* `hdf5`: HDF5 output for checkpoint and visualization files (requires building with `USE_HDF5=1`)
* `adios2`: ADIOS2 output for checkpoint and visualization files (requires building with `USE_ADIOS2=1`)

The `hdf5` and `adios2` backends come from *CaNS*, and have not been tested yet with the two-fluid fields.

```fortran
is_full_restart = F
```

`is_full_restart`, if true, saves in each checkpoint what is needed to **restart the time integration exactly**: the pressure and the right-hand sides of the previous time step (e.g., `fld_po.bin`, `fld_dudtrko.bin`), as one file per field, and the last two time step sizes in `fld_dt.out`. A run restarted from such a checkpoint reproduces an uninterrupted run. By default, or if these files are not found, a restarted run begins with a first-order time step. With the geometric volume-of-fluid method (`INTERFACE_CAPTURING_SPLIT_VOF`), the velocity and the scalar of the previous time step are saved as well (`fld_uo.bin`, `fld_so.bin`, ...).

If `&io` is not provided, *CaNS-Fizzy* defaults to `io_backend = 'mpiio'`. Note also that selecting `hdf5` or `adios2` requires compiling the code with support for the corresponding library; otherwise, *CaNS-Fizzy* falls back to `mpiio` with a warning at runtime.
