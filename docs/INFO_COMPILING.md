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
GPU_BACKEND=OACC   # options: OACC, OMP
USE_HDF5=0         # HDF5 I/O support
USE_ADIOS2=0       # ADIOS2 I/O support
#
# two-fluid options
#
CONSTANT_COEFFS_POISSON=1 # constant-coefficient (FFT-based) pressure solver with pressure splitting
BALANCED_CAPILLARY_PRESSURE_SPLIT=1 # balanced residual of pressure gradient and surface tension in the pressure splitting
```

In this file, `FCOMP` can be one of `GNU` (`gfortran`), `INTEL` (`ifx`), `NVIDIA` (`nvfortran`), or `CRAY` (`ftn`); the predefined profiles for compiler options can be selected by choosing one of the `FFLAGS_*` options; finer control of the compiler flags may be achieved by building with, e.g., `make FFLAGS+=[OTHER_FLAGS]`, or by tweaking the profiles directly under `configs/flags.mk`. Similarly, the library paths (e.g., for *FFTW*) may need to be adapted in the `Makefile` (`LIBS` variable) or by building with `make LIBS+='-L[PATH_TO_LIB] -l[NAME_OF_LIB]'`. Finally, the following pre-processing options are available:

 * `SINGLE_PRECISION` : calculation will be carried out in single precision (the default precision is double); all interface methods can be used, with the geometric one less accurate than the others in this precision (the test `single_precision` prints the differences with respect to double precision for a few cases)
 * `GPU`              : enable GPU accelerated runs
 * `USE_HDF5`         : enable HDF5-based I/O support
 * `USE_ADIOS2`       : enable ADIOS2-based I/O support
 * `USE_DIEZDECOMP`   : use [diezDecomp](https://github.com/Rafael10Diez/diezDecomp) as GPU communication backend instead of the default (cuDecomp), e.g., for portability in different accelerators. While diezDecomp supports CPU-CPU communication, this is not used in CaNS yet.
 * `GPU_BACKEND`      : backend of the directive-based GPU offloading: `OACC` for OpenACC (default), and `OMP` for OpenMP target offload (the OpenMP parallelization on CPUs is then not available)

and the following ones select the two-fluid method (the interface-capturing method itself, ACDI, THINC/QQ or geometric volume of fluid, is selected at run time: see `ifc_reconstruction` in [`INFO_INPUT.md`](INFO_INPUT.md)):

 * `CONSTANT_COEFFS_POISSON`          : uses the direct FFT-based solver for the pressure Poisson equation, with a pressure-splitting method (default); if set to `0`, the variable-coefficient Poisson equation is solved with an iterative multigrid solver of the HYPRE library (CPU builds only; see below)
 * `BALANCED_CAPILLARY_PRESSURE_SPLIT`: extrapolates the balanced residual between the pressure gradient and surface-tension force in the pressure-splitting method, instead of the pressure gradient alone, which jumps across the interface (see P. Cifani, *Journal of Computational Physics* 398 (2019): 108904); it is used by default with `CONSTANT_COEFFS_POISSON=1` (set it to `0` to extrapolate the pressure gradient alone), and has no role with `CONSTANT_COEFFS_POISSON=0`

The transport of a scalar and the Boussinesq buoyancy term are selected in the input file (`nscal` and `is_boussinesq_buoyancy`; see [`INFO_INPUT.md`](INFO_INPUT.md)), and so are the pencil orientation, the debugging checks and the timing of each step.

A build with `CONSTANT_COEFFS_POISSON=0` links to HYPRE, whose location is set with `HYPRE_HOME` (default: `${HOME}/software/hypre`), e.g., `make CONSTANT_COEFFS_POISSON=0 HYPRE_HOME=[PATH_TO_HYPRE]`. The solver is chosen and tuned at run time, with the `&hypre` namelist (see [`INFO_INPUT.md`](INFO_INPUT.md)). Two-dimensional cases should use a single grid cell along the homogeneous direction, which is then left out of the linear system; otherwise, a periodic direction must have an even number of grid cells.

The OpenMP GPU backend is selected by building with `GPU=1 GPU_BACKEND=OMP`; if `GPU_BACKEND` is not specified, the default GPU path remains `OACC`. Each OpenACC directive in the source is followed by its OpenMP counterpart, as in CaNS. With the NVIDIA compiler, the OpenMP build inlines the small procedures called from the kernels (`-Minline`): a loop that calls a procedure is otherwise parallelized across teams only.

NVIDIA GPU builds generate code for the major compute capabilities selected by the compiler (`NVHPC_GPU_TARGET=ccall-major`). This can be overridden for a specific GPU, for example with `make NVHPC_GPU_TARGET=cc89`, or set to `ccnative` to target the GPUs visible on the build system.

Typing `make libs` will build the 2DECOMP&FFT/cuDecomp/diezDecomp libraries; then typing `make` will compile the code and copy the executable `cans` to a `run/` folder; `make run` will also copy the default input file `input.nml` under `src/` to the same `run/` folder.

The executable records the location of the cuDecomp (or hipDecomp) shared libraries, so a GPU run needs no further setup. If the `dependencies` folder is moved after the build, the libraries are found by updating the `LD_LIBRARY_PATH` environment variable instead, e.g. (from the root directory):
```shell
export LD_LIBRARY_PATH=$PWD/dependencies/cuDecomp/build/lib:$LD_LIBRARY_PATH
```

Finally, the choice of compiler `FCOMP` (see `configs/flags.mk`) and profile flags `FFLAGS_*` (see `configs/flags.mk`) can easily be overridden, for instance, as: `make FC=ftn FFLAGS=-O2`. Linking and include options can be changed in `configs/libs.mk`. The default `build.conf` files and `*.mk` files are created from `configs/defaults` at the first compilation.
