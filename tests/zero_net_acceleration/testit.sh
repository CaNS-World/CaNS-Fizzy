#!/bin/bash
#
# needs a build with BALANCED_CAPILLARY_PRESSURE_SPLIT=1
#
set -e # exit immediately if any command fails
#
TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
RUNDIR=$CANSDIR/run

MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi
MPIRUN="mpirun -n 4 $MPIRUN_OPTIONS"

run_case() {
  #
  # runs the case named $1 in its own directory, with the input file $2, and checks it with the test file $3;
  # the optional $4 is a sed expression that edits the input file
  #
  local casedir=$RUNDIR/zero_net_$1
  rm -rf $casedir && mkdir -p $casedir/data
  cp $TESTDIR/$2 $casedir/input.nml && cp $TESTDIR/spheres.in $casedir
  if [ -n "$4" ]; then sed -i "$4" $casedir/input.nml; fi
  echo "INFO: Running CaNS, case $1"
  cd $casedir
  set +e
  ${MPIRUN} $RUNDIR/cans 1> log_file.log 2> err_log.log
  status=$?
  set -e
  cp log_file.log $RUNDIR/zero_net_$1.log
  if [ $status -ne 0 ]; then
    tail -n 50 log_file.log; cat err_log.log
    echo "ERROR: CaNS exited with status $status"
    exit $status
  fi
  cd $casedir/data
  cp $TESTDIR/$3 .
  pytest $3
}

run_case channel input.nml test.py
run_case capillary input_capillary.nml test_capillary.py
run_case capillary_equal_density input_capillary.nml test_capillary_momentum.py 's/rho12(1:2) = 10., 1./rho12(1:2) = 1., 1./'
run_case boussinesq input_boussinesq.nml test_boussinesq.py

#
# inputs that must be rejected
#
check_rejected() {
  #
  # runs the channel case on one task, with the input edited by the sed expressions
  # that follow the case name $1 and the expected error message $2; the run must fail
  #
  local name=$1 message=$2 casedir=$RUNDIR/zero_net_invalid_$1
  shift 2
  rm -rf $casedir && mkdir -p $casedir/data
  cp $TESTDIR/input.nml $casedir/input.nml
  sed -i 's/dims(1:2) = 2, 2/dims(1:2) = 1, 1/' $casedir/input.nml
  for expr in "$@"; do sed -i "$expr" $casedir/input.nml; done
  echo "INFO: Running CaNS, case invalid_$name"
  if (cd $casedir && mpirun -n 1 $RUNDIR/cans > run.log 2>&1); then
    echo "ERROR: the input of case invalid_$name was not rejected"
    exit 1
  fi
  cp $casedir/run.log $RUNDIR/zero_net_invalid_$name.log
  grep -q "$message" $casedir/run.log
}

check_rejected nonperiodic 'requires periodic pressure BCs' \
               "s/'P','P',  'P','P',  'D','D'/'D','D',  'P','P',  'D','D'/" \
               "s/'P','P',  'P','P',  'N','N'/'N','N',  'P','P',  'N','N'/"
check_rejected open 'periodic or impermeable transverse boundaries' \
               's/bcvel(0:1,1:3,3) = 0.,0.,  0.,0.,  0.,0./bcvel(0:1,1:3,3) = 0.,0.,  0.,0.,  1.,1./'
echo "INFO: all cases passed"
