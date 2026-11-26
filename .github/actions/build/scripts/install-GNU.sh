#!/bin/sh
#
# installs the necessary dependencies for building with GNU (`gfortran`) and OpenMPI,
# and HYPRE for the variable-coefficient Poisson solver
#
sudo apt-get update -y -o Acquire::Retries=3
sudo apt-get install -y -o Acquire::Retries=3 gfortran libopenmpi-dev libfftw3-dev libhypre-dev
