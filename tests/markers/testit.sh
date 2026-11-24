#!/bin/bash
#
set -e # exit immediately if any command fails
#
TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
SRCDIR=$CANSDIR/src
RUNDIR=$CANSDIR/run
UTILSDIR=$CANSDIR/utils

MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi
MPIRUN="mpirun -n 4 $MPIRUN_OPTIONS"

cd $RUNDIR
mkdir -p data
#
cp $TESTDIR/input.nml $TESTDIR/spheres-*.in $RUNDIR

#
# with the direction-split advection, the cases are also run with the transport along with the sweeps
#
TRANSPORTS="default"
if [[ "$TEST_TWO_FLUID_INPUT" == *plic* || "$TEST_TWO_FLUID_INPUT" == *split* ]]; then TRANSPORTS="default sweeps"; fi
for transport in $TRANSPORTS; do
  for name in one same split bubbles restart-1 restart-2; do
    echo "INFO: Running CaNS, case $name, transport $transport"
    cp $TESTDIR/input.nml input.nml
    cp spheres-same.in spheres.in
    if [ $name = one   ]; then sed -i "s/^nmarker = 2/nmarker = 1/" input.nml; fi
    if [ $name != one ] && [ $name != same ]; then cp spheres-split.in spheres.in; fi
    #
    # bubbles ten times lighter than the surrounding fluid, instead of drops of the same density
    #
    if [ $name = bubbles ]; then sed -i -e "s/^inipsi.*/inipsi = 'bub2-x'/" -e "s/^rho12.*/rho12(1:2) = 1., 0.1/" input.nml; fi
    #
    # the case `split` in two runs of 300 time steps, restarted exactly
    #
    if [[ $name == restart-* ]]; then
      sed -i -e "s/^nstep = 600/nstep = 300/" -e "s/isave = 600/isave = 300/" input.nml
      printf "\n&io\nis_full_restart = T\n/\n" >> input.nml
    fi
    if [ $name = restart-2 ]; then sed -i -e "s/^nstep = 300/nstep = 600/" -e "s/^restart = F/restart = T/" input.nml; fi
    #
    # interface method of the variant under test, if set (see `tests/testall.sh`): a line for the `two_fluid` namelist
    #
    if [ -n "$TEST_TWO_FLUID_INPUT" ]; then sed -i "/^&two_fluid/a $TEST_TWO_FLUID_INPUT" input.nml; fi
    if [ $transport = sweeps ]; then sed -i "/^&two_fluid/a ifc_transport = 'sweeps'" input.nml; fi
    set +e
    ${MPIRUN} ./cans 1> log_file.log 2> err_log.log
    status=$?
    set -e
    (head -n 50 log_file.log; echo -e "\n[...output omitted...]\n"; tail -n 50 log_file.log) # report first and last N lines of the log file
    if [ -s err_log.log ]; then
      echo -e "\nINFO: stderr from CaNS\n"
      cat err_log.log
    fi
    if [ $status -ne 0 ]; then
      echo "ERROR: CaNS exited with status $status"
      exit $status
    fi
    if [ $name = restart-1 ]; then continue; fi # (its checkpoint stays in place for the second run)
    mkdir -p data-$name
    mv data/* data-$name
  done

  cp $TESTDIR/*.* .
  echo "INFO: Checking the results"
  pytest test.py
  rm -rf $RUNDIR/data/*.* data-*
done
