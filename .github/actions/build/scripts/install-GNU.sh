#!/bin/sh
#
# installs the necessary dependencies for building with GNU (`gfortran`) and OpenMPI,
# and HYPRE for the variable-coefficient Poisson solver
#
sudo apt-get install gfortran libopenmpi-dev libfftw3-dev libhypre-dev
