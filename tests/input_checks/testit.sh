#!/bin/bash
#
# the default input file of the code, `src/input.nml` (the one that `make run` copies), must run as it is,
# and inputs with mistakes that the code checks must be rejected with the corresponding message
#
set -e # exit immediately if any command fails
#
TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
SRCDIR=$CANSDIR/src
RUNDIR=$CANSDIR/run

MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi
MPIRUN="mpirun -n 4 $MPIRUN_OPTIONS"

cd $RUNDIR
mkdir -p data
cp $SRCDIR/input.nml $SRCDIR/spheres.in $RUNDIR
echo "INFO: Running CaNS, default input file"
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
grep -q "Fim" log_file.log
rm -rf $RUNDIR/data/*.*

#
# inputs that must be rejected
#
check_rejected() {
  #
  # runs the default input on one task, edited with the sed expression $3, and with $4 (if any) as the
  # contents of `spheres.in`; the run must fail with the message $2
  #
  local name=$1 message=$2 expr=$3 spheres=$4 casedir=$RUNDIR/input_checks_$1
  rm -rf $casedir && mkdir -p $casedir/data
  sed -e "s/^dims(1:2) *=.*/dims(1:2) = 1, 1/" -e "$expr" $SRCDIR/input.nml > $casedir/input.nml
  if [ -n "$spheres" ]; then printf "$spheres" > $casedir/spheres.in; fi
  echo "INFO: Running CaNS, case invalid_$name"
  if (cd $casedir && mpirun -n 1 $RUNDIR/cans > run.log 2>&1); then
    echo "ERROR: the input of case invalid_$name was not rejected"
    exit 1
  fi
  cp $casedir/run.log $RUNDIR/input_checks_$name.log
  grep -q "$message" $casedir/run.log
  rm -rf $casedir
}

check_rejected psi_periodicity 'must match pressure BCs' "s/^cbcpsi.*/cbcpsi(0:1,1:3) = 'P','P',  'P','P',  'N','N'/"
check_rejected psi_pair        'not valid'               "s/^cbcpsi.*/cbcpsi(0:1,1:3) = 'P','N',  'P','P',  'P','P'/"
check_rejected density         'must be positive'        "s/^rho12.*/rho12(1:2) = 1., 0./"
check_rejected marker_range    'must be between 1 and'   "/^&two_fluid/a nmarker = 2" "0.5 0.5 0.3 0.15 1\n0.5 0.5 0.7 0.15 3\n"
check_rejected marker_file     'need the file'           "/^&two_fluid/a nmarker = 2"
echo "INFO: all cases passed"
