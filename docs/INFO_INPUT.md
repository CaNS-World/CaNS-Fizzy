# about the input file `input.nml`

Consider the following input file as example (corresponds to a turbulent plane channel flow with bubbles and walls at different temperatures). `&dns` defines a so-called Fortran namelist containing the necessary computational parameters to set a case.


```fortran
&dns
ng(1:3) = 512, 256, 144
l(1:3) = 6., 3., 1.
gtype = 1, gr = 0.
cfl = 0.95, dtmin = 1.e5, dt_f = -1
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
gacc(1:3) = 0., 0., 0.
dims(1:2) = 2, 2
\
```

Additionally, the `&scalar` namelist contains the computational parameters for scalar transport (e.g. temperature), while the physical and transport properties of the two phases and the numerical parameters of the two-phase solver are specified in the `&two_fluid` namelist.

```fortran
&scalar
inisca          = 'dht'
cbcsca(0:1,1:3) = 'P','P',  'P','P',  'D','D'
bcsca(0:1,1:3)  =  0.,0. ,   0.,0. ,   1.,0.
ssource         = 0.
\

&two_fluid
inipsi            = 'bub3'
cbcpsi(0:1,1:3)   = 'P','P',  'P','P',  'P','P'
cbcnor(0:1,1:3,1) = 'P','P',  'P','P',  'P','P'
cbcnor(0:1,1:3,2) = 'P','P',  'P','P',  'P','P'
cbcnor(0:1,1:3,3) = 'P','P',  'P','P',  'P','P'
bcpsi(0:1,1:3)    =  0.,0. ,   0.,0. ,   0.,0.
bcnor(0:1,1:3,1)  =  0.,0.,    0.,0.,    0.,0.
bcnor(0:1,1:3,2)  =  0.,0.,    0.,0.,    0.,0.
bcnor(0:1,1:3,3)  =  0.,0.,    0.,0.,    0.,0.
sigma             =  1.
rho12(1:2)        =  1., 1.
mu12(1:2)         =  100., 100.
ka12(1:2)         =  0., 0.
cp12(1:2)         =  0., 0.
beta12(1:2)       =  0., 0.
acdi_gam_factor   =  1., acdi_gam_min = 0.
acdi_eps_factor   =  1.
\
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

`gtype` and `gr` are the **grid stretching type** and **grid stretching parameter** that tweak the non-uniform grid in the third direction; zero `gr` implies no stretching. See `initgrid.f90` for more details. The following options are available for `gtype`:

* `1`: grid clustered towards both ends (default)
* `2`: grid clustered towards the lower end
* `3`: grid clustered towards the upper end
* `4`: grid clustered towards the middle

---

```fortran
cfl = 0.95, dtmin = 1.e5, dt_f = -1
```

This line controls the simulation time step.

The time step is set to be equal to `min(cfl*dtmax,dtmin)`, i.e. the minimum value between `dtmin` and `cfl` times the maximum allowable time step `dtmax` (computed every `ickeck` time steps; see below).
`dtmin` is therefore used when a constant time step, smaller than `cfl*dtmax`, is required. If not, it should be set to a high value so that the time step is dynamically adjusted to `cfl*dtmax`. Alternatively, one can force the simulation to advance with a constant time step of arbitrary value, independent of the `dtmax` evaluation, by changing the value of `dt_f` from `-1` to the desired time step.

---

```fortran
is_solve_ns = T, is_track_interface = T
```

`is_solve_ns`, if true, enables the solution of the Navier-Stokes momentum equation. It can be set to false when a prescribed velocity field is desired.
`is_track_interface`, if true, enables the two-phase interface advection.

---

```fortran
inivel = 'poi'
is_wallturb = T, is_forced_hit = F
```

These lines set the initial velocity field.

`initvel` **chooses the initial velocity field**. The following options are available:

* `zer`: zero velocity field
* `uni`: uniform velocity field equal to `uref`                                     ; streamwise direction in `x`
* `cou`: plane Couette flow profile with symmetric wall velocities equal to `uref/2`; streamwise direction in `x`
* `poi`: plane Poiseuille flow profile with mean velocity `uref`                    ; streamwise direction in `x`
* `iop`: same as `poi`, with the reference frame moving with velocity `uref`        ; streamwise direction in `x`
* `tbl`: temporal boundary layer profile with wall velocity `uref`                  ; streamwise direction in `x`
* `pdc`: plane Poiseuille flow profile with constant pressure gradient              ; streamwise direction in `x`
* `log`: logarithmic profile with mean velocity `uref`                              ; streamwise direction in `x`
* `hcp`: half channel with plane Poiseuille profile and mean velocity `uref`        ; streamwise direction in `x`
* `hcl`: half channel with logarithmic profile and mean velocity `uref`             ; streamwise direction in `x`
* `hdc`: half plane Poiseuille flow profile with constant pressure gradient         ; streamwise direction in `x`
* `tgv`: three-dimensional Taylor-Green vortex
* `tgw`: two-dimensional   Taylor-Green vortex
* `rox`: solid body rotation with angular velocity $2 \pi$                          ; around the `x` axis
* `roy`: solid body rotation with angular velocity $2 \pi$                          ; around the `y` axis
* `roz`: solid body rotation with angular velocity $2 \pi$                          ; around the `z` axis
* `vox`: two-dimensional, time dependent serpentine vortex                          ; around the `x` axis
* `voy`: two-dimensional, time dependent serpentine vortex                          ; around the `y` axis
* `voz`: two-dimensional, time dependent serpentine vortex                          ; around the `z` axis

`is_wallturb`, if true, **superimposes a high amplitude disturbance on the initial velocity field** that effectively triggers transition to turbulence in a wall-bounded shear flow.

See `initflow.f90` for more details.

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

a checkpoint file `fld.bin` will be saved before the simulation is terminated.

`restart`, if true, **restarts the simulation** from a previously saved checkpoint file, named `fld.bin`.

`is_overwrite_save`, if true, overwrites the checkpoint file `fld.bin` at every save; if false, a symbolic link is created which makes `fld.bin` point to the last checkpoint file with name `fld_???????.bin` (with `???????` denoting the corresponding time step number). In the latter case, to restart a run from a different checkpoint one just has to point the file `fld.bin` to the right file, e.g.: ` ln -sf fld_0000100.bin fld.bin`.

`nsaves_max` limits the number of saved checkpoint files, if `is_over_write_save` is false; a value of `0` or any negative integer corresponds to no limit, and the code uses the file format described above; otherwise, only `nsaves_max` checkpoint files are saved, with the oldest save being overwritten when the number of saved checkpoints exceeds this threshold; in this case, files with a format `fld_????.bin` are saved (with `????` denoting the saved file number), with `fld.bin` pointing to the last checkpoint file as described above; moreover, a file `log_checkpoints.out` records information about the time step number and physical time corresponding to each saved file number.

---

```fortran
icheck = 10, iout0d = 10, iout1d = 100, iout2d = 500, iout3d = 10000, isave = 5000
```

These lines set the frequency of time step checking and output:

* every `icheck` time steps **the new time step size** is computed according to the new stability criterion and cfl (above)
* every `iout0d` time steps **history files with global scalar variables** are appended
* every `iout1d` time steps **1d profiles** are written (e.g. velocity and its moments) to a file
* every `iout2d` time steps **2d slices of a 3d scalar field** are written to a file
* every `iout3d` time steps **3d scalar fields** are written to a file
* every `isave`  time steps a **checkpoint file** is written (`fld_???????.bin`), and a symbolic link for the restart file, `fld.bin`, will point to this last save so that, by default, the last saved checkpoint file is used to restart the simulation

1d, 2d and 3d outputs can be tweaked modifying files `out?d.h90`, and re-compiling the source. See also `output.f90` for more details.

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

These lines set the boundary conditions (BC) for velocity and pressure.

The **type** (BC) for each field variable is set by a row of six characters, `X0 X1  Y0 Y1  Z0 Z1` where,

* `X0` `X1` set the type of BC for the field variable on the **lower** and **upper** boundaries in `x`
* `Y0` `Y1` set the type of BC for the field variable on the **lower** and **upper** boundaries in `y`
* `Z0` `Z1` set the type of BC for the field variable on the **lower** and **upper** boundaries in `z`

The four rows correspond to the three velocity components, and pressure, i.e. `u`, `v`, `w`, and `p`.

The following options are available:

* `P` periodic
* `D` Dirichlet
* `N` Neumann

The **last four rows** follow the same logic, but now for the BC **values** (dummy for a periodic direction).

---

```fortran
bforce(1:3) = 0., 0., 0.
gacc(1:3) = 0., 0., 0.
```

These lines set the flow forcing.

`bforce`, is a constant **body force density term** in the direction in question (e.g., the negative of a constant pressure gradient) that can be added to the right-hand-side of the momentum equation. The three values correspond to three domain directions.
`gacc`, is a constant **acceleration** in the direction in question (e.g., gravity) added to the right-hand-side of the momentum equation. The three values correspond to three domain directions.

---

```fortran
dims(1:2) = 2, 2
```

This line sets the grid of computational subdomains.

`dims` is the **processor grid**, the number of domain partitions along the first and second decomposed directions (which depend on the selected default pencil orientation). `dims(1)*dims(2)` corresponds therefore to the total number of computational subdomains. Setting `dims(:) = [0,0]` will trigger a runtime autotuning step to find the processor grid that minimizes transpose times. Note, however, that other components of the algorithm (e.g., collective I/O) may also be affected by the choice of processor grid.

---

```fortran
inisca          = 'dht'
```

`inisca` **chooses the initial scalar field** (e.g. temperature). The following options are available:

* `zer`: zero scalar field
* `uni`: uniform scalar field equal to `sref`
* `cou`: plane scalar profile, linear between the wall boundary values
* `dhc`: differentially heated cavity, in the x direction

See `initflow.f90` for more details.

---

```fortran
cbcsca(0:1,1:3) = 'P','P',  'P','P',  'P','P'
bcsca(0:1,1:3)  =  0.,0. ,   0.,0. ,   0.,0.
```

These lines set the boundary conditions (BC) for the scalar.

The **type** (BC) is set in `cbcsca` by a row of six characters, `X0 X1  Y0 Y1  Z0 Z1` where,

* `X0` `X1` set the type of BC for the scalar on the **lower** and **upper** boundaries in `x`
* `Y0` `Y1` set the type of BC for the scalar on the **lower** and **upper** boundaries in `y`
* `Z0` `Z1` set the type of BC for the scalar on the **lower** and **upper** boundaries in `z`

The following options are available:

* `P` periodic
* `D` Dirichlet
* `N` Neumann

The corresponding BC **values** are set in `bcsca` (dummy for a periodic direction).

---

```fortran
ssource         = 0.
```

This line sets the value of a scalar source, that can be added to the right-hand-side of the scalar transport equation.

---

```fortran
inipsi          = 'bub3'
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
* `cap-wav-1d`: planar small-amplitude capillary wave
* `zalesak-disk`: two-dimensional (lighter phase) slotted disk

See `two_fluid.f90` for more details. For the `bub` and `drp` initial fields, the position and size of the films/bubbles/droplets are specified in the `spheres.in` file, where each line corresponds to an individual bubble/droplet, specifying the cartesian coordinates of its center of mass, and its radius.
For example, the following file introduces in the computational domain two bubbles/droplets, the first centered at `x=0.3`, `y=1.`, `z=0.5` and with a radius of `0.25`, the second centered at `x=1.5`, `y=0.6`, `z=1.` and with a radius of `0.5`.

```fortran
0.3 1.  0.5 0.25
1.5 0.6 1.  0.5
```

If `spheres.in` is used to initialize a planar film (options `bub1` and `drp1`), then the radius value will represent the film thickness, and the film will be aligned along the XY plane.

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
For consistency, the components of the interface normal should be periodic when the volume fraction boundary conditions are periodic. In order to conserve mass in the domain, a wall boundary requires a homogenous Neumann BC for the volume fraction and a homogeneous Dirichlet BC for the wall-normal component of the interface normal. In this case the BC for the wall-tangential components of the interface normal should be set to homogenous Neumann.

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

`sigma` sets the value of the surface tension at the two-phase interface.
`rho12` sets the value of the density of each phase. _Note: the initial conditions have been designed to consider phase 1 as the densest one._
`mu12` sets the value of the dynamic viscosity of each phase.
`ka12` sets the value of the thermal conductivity of each phase.
`cp12` sets the value of the constant pressure specfic heat capacity of each phase.
`beta12` sets the value of the thermal expansion coeffcient of each phase.

---

```fortran
acdi_gam_factor = 1., acdi_gam_min = 0.
acdi_eps_factor = 0.51
```

These lines set the computational parameters specific to the phase field solver.

`acdi_gam_factor` sets the speed of the interface relaxation, relative to the maximum flow velocity. It should always be equal to or larger than `1`.
`acdi_gam_min` sets a lower threshold for the speed of the interface relaxation.
`acdi_eps_factor` sets the interface thickness relative to the maximum grid size. It should always be larger than `0.5`.


# about the `&cudecomp` namelist under `input.nml`

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
