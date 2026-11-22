#!/bin/bash
#
# builds and runs the unit tests of the kernels of the PLIC and THINC reconstructions (`src/plic.f90` and
# `src/thinc.f90`), which do not need the solver;
# the compiler and its flags can be set with the environment variables FC and FFLAGS (default: GNU)
#
set -e # exit immediately if any command fails
#
TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
SRCDIR=$CANSDIR/src
RUNDIR=$CANSDIR/run
FC=${FC:-mpifort}
FFLAGS=${FFLAGS:--O2 -cpp -Wall -Wextra -fcheck=all}

BUILDDIR=$RUNDIR/plic_unit
rm -rf $BUILDDIR && mkdir -p $BUILDDIR && cd $BUILDDIR
echo "INFO: Building the unit tests"
$FC $FFLAGS -c $SRCDIR/types.f90 $SRCDIR/param.f90 $SRCDIR/plic.f90 $SRCDIR/thinc.f90
$FC $FFLAGS $TESTDIR/test_plic.f90 types.o param.o plic.o thinc.o -o test_plic
echo "INFO: Running the unit tests"
./test_plic
python3 $TESTDIR/check_exact.py volume_samples.txt
cd $RUNDIR && rm -rf $BUILDDIR
