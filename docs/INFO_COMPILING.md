# Compiling CaNS-Fizzy

For most systems, CaNS-Fizzy can be compiled from the root directory with the following commands `make libs && make`, which will compile the 2DECOMP&FFT/cuDecomp libraries and CaNS-Fizzy. `make clean` clears the CaNS-Fizzy build files, `make libsclean` clears the 2DECOMP/cuDecomp builds, and `make allclean` clears both.

The `Makefile` in the root directory is used to compile the code, and is expected to work out-of-the-box for most systems. The `build.conf` file in the root directory can be used to choose the Fortran compiler (MPI wrapper), a few pre-defined profiles depending on the nature of the run (e.g., production vs debugging), and pre-processing options:

```shell
#
# compiler and compiling profile
#
FCOMP=GNU          # options: GNU, NVIDIA, INTEL
FFLAGS_OPT=1       # for production runs
FFLAGS_OPT_MAX=0   # for production runs (more aggressive optimization)
FFLAGS_DEBUG=0     # for debugging
FFLAGS_DEBUG_MAX=0 # for thorough debugging
#
# defines
#
SINGLE_PRECISION=0 # perform the whole calculation in single precision
GPU=0              # GPU build
USE_HDF5=0         # HDF5 I/O support
USE_ADIOS2=0       # ADIOS2 I/O support
#
# two-fluid options
#
CONSTANT_COEFFS_POISSON=1  # constant-coefficient (FFT-based) pressure solver with pressure splitting
INTERFACE_CAPTURING_VOF=0  # uses THINC/QQ volume-of-fluid instead of the ACDI pathway
INTERFACE_CAPTURING_SPLIT_VOF=0 # uses the geometric (PLIC) direction-split volume-of-fluid instead of the ACDI pathway
SDF_NORMALS=1              # default: ACDI or SDF-based THINC reconstruction and normals
```

In this file, `FCOMP` can be one of `GNU` (`gfortran`), `INTEL` (`ifx`), `NVIDIA` (`nvfortran`), or `CRAY` (`ftn`); the predefined profiles for compiler options can be selected by choosing one of the `FFLAGS_*` options; finer control of the compiler flags may be achieved by building with, e.g., `make FFLAGS+=[OTHER_FLAGS]`, or by tweaking the profiles directly under `configs/flags.mk`. Similarly, the library paths (e.g., for *FFTW*) may need to be adapted in the `Makefile` (`LIBS` variable) or by building with `make LIBS+='-L[PATH_TO_LIB] -l[NAME_OF_LIB]'`. Finally, the following pre-processing options are available:

 * `SINGLE_PRECISION` : calculation will be carried out in single precision (the default precision is double)
 * `GPU`              : enable GPU accelerated runs
 * `USE_HDF5`         : enable HDF5-based I/O support
 * `USE_ADIOS2`       : enable ADIOS2-based I/O support
 * `USE_DIEZDECOMP`   : use [diezDecomp](https://github.com/Rafael10Diez/diezDecomp) as GPU communication backend instead of the default (cuDecomp), e.g., for portability in different accelerators. While diezDecomp supports CPU-CPU communication, this is not used in CaNS yet.

and the following ones select the two-fluid method:

 * `CONSTANT_COEFFS_POISSON`          : uses the direct FFT-based solver for the pressure Poisson equation, with a pressure-splitting method (default); if set to `0`, the variable-coefficient Poisson equation is solved with an iterative multigrid solver of the HYPRE library (CPU builds only; see below)
 * `INTERFACE_CAPTURING_VOF`          : uses the THINC/QQ volume-of-fluid method instead of the ACDI diffuse-interface pathway
 * `INTERFACE_CAPTURING_SPLIT_VOF`    : uses the geometric volume-of-fluid method instead of the ACDI diffuse-interface pathway: a piecewise-linear (PLIC) reconstruction of the interface, and a conservative direction-split advection that conserves the volume and keeps the volume fraction bounded to round-off, with sub-steps when the CFL number along a direction exceeds 1/2. For now, this method only transports the interface in a prescribed velocity field (`is_solve_ns = F`, without a scalar); it cannot be combined with `INTERFACE_CAPTURING_VOF`, and it does not use `SDF_NORMALS`
 * `SDF_NORMALS`                      : enables approximate signed-distance fields for interface transport, normals, and curvature (default); setting it to `0` uses CDI or standard THINC directly from the phase indicator without allocating a distance field
 * `BALANCED_CAPILLARY_PRESSURE_SPLIT`: extrapolates the balanced residual between the pressure gradient and surface-tension force in the pressure-splitting method
 * `CAPILLARY_BRACKBILL_NORMALIZATION`: enables density-based Brackbill scaling of the surface-tension force (default); setting it to `0` disables the scaling
 * `VISC_HARMONIC_INTERPOLATION`      : uses harmonic phase interpolation for viscosity and thermal conductivity instead of the default arithmetic interpolation

The transport of a scalar and the Boussinesq buoyancy term are selected in the input file (`nscal` and `is_boussinesq_buoyancy`; see [`INFO_INPUT.md`](INFO_INPUT.md)), and so are the pencil orientation, the debugging checks and the timing of each step.

A build with `CONSTANT_COEFFS_POISSON=0` links to HYPRE, whose location is set with `HYPRE_HOME` (default: `${HOME}/software/hypre`), e.g., `make CONSTANT_COEFFS_POISSON=0 HYPRE_HOME=[PATH_TO_HYPRE]`. The solver is chosen and tuned at run time, with the `&hypre` namelist (see [`INFO_INPUT.md`](INFO_INPUT.md)). Two-dimensional cases should use a single grid cell along the homogeneous direction, which is then left out of the linear system; otherwise, a periodic direction must have an even number of grid cells.

The GPU-enabled path corresponds to the OpenACC implementation.

NVIDIA GPU builds generate code for the major compute capabilities selected by the compiler (`NVHPC_GPU_TARGET=ccall-major`). This can be overridden for a specific GPU, for example with `make NVHPC_GPU_TARGET=cc89`, or set to `ccnative` to target the GPUs visible on the build system.

Typing `make libs` will build the 2DECOMP&FFT/cuDecomp/diezDecomp libraries; then typing `make` will compile the code and copy the executable `cans` to a `run/` folder; `make run` will also copy the default input file `input.nml` under `src/` to the same `run/` folder.

The executable records the location of the cuDecomp (or hipDecomp) shared libraries, so a GPU run needs no further setup. If the `dependencies` folder is moved after the build, the libraries are found by updating the `LD_LIBRARY_PATH` environment variable instead, e.g. (from the root directory):
```shell
export LD_LIBRARY_PATH=$PWD/dependencies/cuDecomp/build/lib:$LD_LIBRARY_PATH
```

Finally, the choice of compiler `FCOMP` (see `configs/flags.mk`) and profile flags `FFLAGS_*` (see `configs/flags.mk`) can easily be overridden, for instance, as: `make FC=ftn FFLAGS=-O2`. Linking and include options can be changed in `configs/libs.mk`. The default `build.conf` files and `*.mk` files are created from `configs/defaults` at the first compilation.
