#!/bin/bash
#
# a single fluid without interface tracking (`is_track_interface = F`) must give the same flow whether
# it is phase 1 (`inipsi = 'uni'`) or phase 2 (`'zer'`), when the two phases have different properties
#
set -e # exit immediately if any command fails
#
TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
RUNDIR=$CANSDIR/run

MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi
MPIRUN="mpirun -n 4 $MPIRUN_OPTIONS"

cd $RUNDIR
mkdir -p data
#
# the fluid has a density of 2 and a viscosity of 0.02, as phase 1 or as phase 2, with each transport of
# momentum: AB2, and with the direction-split advection, in one step and along with the sweeps
#
for transport in ab2 one-step sweeps; do
  for phase in 1 2; do
    echo "INFO: Running CaNS, transport $transport, phase $phase"
    cp $TESTDIR/input.nml input.nml
    if [ $phase = 2 ]; then
      sed -i -e "s/^inipsi.*/inipsi = 'zer'/" -e "s/^rho12.*/rho12(1:2) = 1., 2./" -e "s/^mu12.*/mu12(1:2) = 0.01, 0.02/" input.nml
    fi
    if [ $transport != ab2 ]; then
      sed -i "/^&two_fluid/a ifc_reconstruction = 'plic', ifc_transport = '$transport'" input.nml
    fi
    set +e
    ${MPIRUN} ./cans 1> log_file.log 2> err_log.log
    status=$?
    set -e
    if [ -s err_log.log ]; then
      echo -e "\nINFO: stderr from CaNS\n"
      cat err_log.log
    fi
    if [ $status -ne 0 ]; then
      tail -n 50 log_file.log
      echo "ERROR: CaNS exited with status $status"
      exit $status
    fi
    mkdir -p data-$transport-$phase
    mv data/* data-$transport-$phase
  done
done

cp $TESTDIR/*.* .
echo "INFO: Checking the results"
pytest test.py
rm -rf $RUNDIR/data/*.* data-*
