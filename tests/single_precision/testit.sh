#!/bin/bash
#
# compares single- with double-precision runs of a few short cases: the single-precision executable is built
# here, in a copy of the source tree (the build of the tree, in double precision, is left as it is), and the
# cases are those of other tests, each run with the two executables
#
set -e # exit immediately if any command fails
#
TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
TESTSDIR=$CANSDIR/tests
RUNDIR=$CANSDIR/run
BUILDDIR=$RUNDIR/single_precision

MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi
MPIRUN="mpirun -n 4 $MPIRUN_OPTIONS"
#
# the comparison is meant for optimized builds: skipped in the short mode of the slow builds (see `testall.sh`)
#
if [ -n "$TEST_SHORT" ]; then echo "INFO: single-precision comparison skipped in the short test mode"; exit 0; fi

echo "INFO: Building CaNS in single precision"
rm -rf $BUILDDIR && mkdir -p $BUILDDIR/dependencies
cd $CANSDIR
cp -r Makefile configs src $BUILDDIR
if [ -f build.conf ]; then cp build.conf $BUILDDIR; fi
cp -r dependencies/2decomp-fft dependencies/configs dependencies/external.mk $BUILDDIR/dependencies
cd $BUILDDIR
rm -rf build run
set +e
(make libsclean && make libs SINGLE_PRECISION=1 && make -j ${NJOBS:-4} SINGLE_PRECISION=1) > $RUNDIR/build_single.log 2>&1
status=$?
set -e
if [ $status -ne 0 ]; then
  tail -n 50 $RUNDIR/build_single.log
  echo "ERROR: the single-precision build failed"
  exit $status
fi
#
# cases: name, input file, file with the spheres (- if none), number of time steps,
#        and line for the `two_fluid` namelist (- if none)
#
CASES=(
"bubble-acdi       rising_bubble/input.nml      rising_bubble/spheres.in               2000 ifc_reconstruction='acdi'"
"bubble-thinc      rising_bubble/input.nml      rising_bubble/spheres.in               2000 ifc_reconstruction='thinc'"
"bubble-plic       rising_bubble/input.nml      rising_bubble/spheres.in               2000 ifc_reconstruction='plic'"
"bubble-thincsplit rising_bubble/input.nml      rising_bubble/spheres.in               2000 ifc_reconstruction='thinc',ifc_advection='split'"
"drop-plic         contact_angle/input-plic.nml contact_angle/spheres.in               2000 -"
"advection-plic    plic_advection/input-3d.nml  plic_advection/spheres-3d.in            768 -"
"markers-plic      markers/input.nml            markers/spheres-split.in                600 ifc_reconstruction='plic',ifc_transport='sweeps'"
"marangoni-plic    marangoni_markers/input.nml  marangoni_markers/spheres-pair-split.in 100 -"
)
cd $RUNDIR
mkdir -p data
for case in "${CASES[@]}"; do
  read -r name input spheres nstep line <<< "$case"
  for precision in double single; do
    echo "INFO: Running CaNS, case $name, $precision precision"
    #
    # (a fixed number of time steps, and a checkpoint at the end)
    #
    sed -e "s/^nstep *= *[0-9]*/nstep = $nstep/" -e "s/^stop_type(1:3) *=.*/stop_type(1:3) = T, F, F/" \
        -e "s/isave *= *[0-9]*/isave = $nstep/" $TESTSDIR/$input > input.nml
    if [ "$line" != "-" ]; then sed -i "/^&two_fluid/a $line" input.nml; fi
    if [ "$spheres" != "-" ]; then cp $TESTSDIR/$spheres spheres.in; fi
    exe=./cans
    if [ $precision = single ]; then exe=$BUILDDIR/run/cans; fi
    set +e
    ${MPIRUN} $exe 1> log_file.log 2> err_log.log
    status=$?
    set -e
    tail -n 20 log_file.log
    if [ -s err_log.log ]; then
      echo -e "\nINFO: stderr from CaNS\n"
      cat err_log.log
    fi
    if [ $status -ne 0 ]; then
      echo "ERROR: CaNS exited with status $status"
      exit $status
    fi
    mkdir -p data-$name-$precision
    mv data/* data-$name-$precision
  done
done

cp $TESTDIR/*.* .
echo "INFO: Checking the results"
pytest -s test.py
rm -rf $RUNDIR/data/*.* data-* $BUILDDIR
