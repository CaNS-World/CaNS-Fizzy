## Synopsis

**CaNS-Fizzy** is a code for massively-parallel numerical simulations of two-phase flows. It has been designed to be an efficient two-fluid Navier-Stokes solver taking [CaNS](https://github.com/CaNS-World/CaNS) as its base, where it is ensured that important advancements in the base solver are easily incorporated. The code aims at solving any fluid flow of two immiscible, incompressible, Newtonian fluid phases that can benefit from a FFT-based solver for the second-order finite-difference Poisson equation in a 3D Cartesian grid. The two-fluid Navier-Stokes equations are solved using a pressure-splitting technique that enables the usage of fast solvers for constant-coefficients Poisson equations. The interface between phases is captured using an Accurate Conservative Diffuse Interface (ACDI) method. See the references below for more details.

**References**

P. Costa. *A FFT-based finite-difference solver for massively-parallel direct numerical simulations of turbulent flows.* *Computers & Mathematics with Applications* 76: 1853--1862 (2018). [doi:10.1016/j.camwa.2018.07.034](https://doi.org/10.1016/j.camwa.2018.07.034) [[arXiv preprint]](https://arxiv.org/abs/1802.10323)

G. Frantzis, & D. Grigoriadis. *An efficient method for two-fluid incompressible flows appropriate for the immersed boundary method.* *Journal of Computational Physics* 376 (2019): 28-53. [doi.org/10.1016/j.jcp.2018.09.035](https://doi.org/10.1016/j.jcp.2018.09.035).

S. Jain. *Accurate conservative phase-field method for simulation of two-phase flows.* *Journal of Computational Physics* 469: 111529 (2022). [doi.org/10.1016/j.jcp.2022.111529](https://doi.org/10.1016/j.jcp.2022.111529)

## Features

Some features are:

 * MPI parallelization
 * FFTW guru interface / cuFFT used for computing multi-dimensional vectors of 1D transforms
 * The right type of transformation (Fourier, cosine, sine, etc) is automatically determined form the input file
 * [cuDecomp](https://github.com/NVIDIA/cuDecomp) pencil decomposition library for _hardware-adaptive_ distributed memory calculations on _many GPUs_
 * [2DECOMP&FFT](https://github.com/xcompact3d/2decomp-fft) library used for performing global data transpositions on CPUs and some of the data I/O
 * GPU acceleration using OpenACC directives
 * A different canonical flow can be simulated just by changing the input files

## Motivation

This numerical tool serves as an efficient base two-fluid Navier-Stokes solver for high-resolution simulations of turbulent multiphase flows. It enables the simulation of two-phase flow in canonical configurations on modern (GPU-based) parallel computing architectures, taking advantage of the efficiency of the base single-phase solver that was used as the starting point. The ACDI method for interface capturing has been selected due to its accuracy and suitability for GPU-based simulations. The one-fluid formulation, combined with the one-equation transport of the phase field, facilitates the extension of the code to more complex numerical strategies to, for instance, introduce complex geometries (such as immersed-boundary methods).

## Method

The two-phase flow is described by a one-fluid formulation, and is solved with a second-order finite difference incremental pressure correction scheme, discretized in a staggered grid arrangement. The interface between the two fluid phases is represented by a diffuse interface of specified thickness, preserved by a regularization flux, and advected with a second-order finite difference scheme. Time is advanced with a three-step low storage Runge-Kutta scheme, with several possible options concerning the discretiation of the advective terms. A pressure splitting technique is used to convert the problem of solving a variable-coefficients Poisson equation to a constant-coefficients one that can leverage the fast Poisson solver in CaNS. See the references above for details.

## Usage

### Downloading *CaNS-Fizzy*

Since *CaNS-Fizzy* loads the external pencil decomposition libraries as Git Submodules, the repository should be cloned as follows:
```bash
git clone --recursive https://github.com/CaNS-World/CaNS-Fizzy
```
so the libraries are downloaded too. Alternatively, in case the repository has already been cloned without the Submodules (i.e., folders `cuDecomp` and `2decomp-fft` under `dependencies/` are empty), the following command can be used to update them:
```bash
git submodule update --init --recursive
```

### Compilation

#### Prerequisites
The prerequisites for compiling CaNS-Fizzy are the following:

 * MPI
 * FFTW3/cuFFT library for CPU/GPU runs
 * The `nvfortran` compiler (for GPU runs)
 * NCCL and NVSHMEM (optional, may be exploited by the cuDecomp library)
 * HYPRE library in case the variable-coefficients Poisson equation is solved without the pressure splitting technique (only available for CPU runs)

#### In short
For most systems, CaNS-Fizzy can be compiled from the root directory with the following commands `make libs && make`, which will compile the 2DECOMP&FFT/cuDecomp libraries, and CaNS-Fizzy.

#### Detailed instructions
The `Makefile` in root directory is used to compile the code, and is expected to work out-of-the-box for most systems. The `build.conf` file in the root directory can be used to choose the Fortran compiler (MPI wrapper), and a few pre-defined profiles depending on the nature of the run (e.g., production vs debugging), and pre-processing options. The following general pre-processing options are available:

 * `DEBUG`                    : performs some basic checks for debugging purposes
 * `TIMING`                   : wall-clock time per time step is computed
 * `PENCIL_AXIS`              : sets the default pencil direction, one of [1,2,3] for [X,Y,Z]-aligned pencils; X-aligned is the default and should be optimal for all cases
 * `SINGLE_PRECISION`         : calculation will be carried out in single precision (the default precision is double)
 * `GPU`                      : enable GPU-accelerated runs

See [`INFO_COMPILING.md`](docs/INFO_COMPILING.md) for more compilation details and a comprehensive description of all pre-processing options.

### Input file

The input file `input.nml` sets the physical and computational parameters. In the `examples/` folder are examples of input files for several canonical flows. See [`INFO_INPUT.md`](docs/INFO_INPUT.md) for a detailed description of the input file.

Files `out1d.h90`, `out2d.h90` and `out3d.h90` in `src/` set which data are written in 1-, 2- and 3-dimensional output files, respectively. *The code should be recompiled after editing out?d.h90 files*.

### Running the code

Run the executable with `mpirun` with a number of tasks complying to what has been set in the input file `input.nml`. Data will be written by default in a folder named `data/`, which must be located where the executable is run (by default in the `run/` folder).

### Visualizing field data

See [`INFO_VISU.md`](docs/INFO_VISU.md).

## Contributing

We appreciate any contributions and feedback that can improve the code. If you wish to contribute to the tool, please get in touch with the maintainers or open an Issue in the repository / a thread in Discussions. Pull Requests are welcome, but please propose/discuss the changes in an linked Issue first.

## Final notes

Please read the `ACKNOWLEDGEMENTS`, `LICENSE` files.
