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
cp $TESTDIR/input-*.nml $TESTDIR/spheres.in $RUNDIR
#
# interface method of the variant under test, if set (see `tests/testall.sh`): a line for the `two_fluid` namelist
#
if [ -n "$TEST_TWO_FLUID_INPUT" ]; then sed -i "/^&two_fluid/a $TEST_TWO_FLUID_INPUT" input*.nml; fi

NAME_EXP="explicit"
NAME_IMP_Z="implicit-z"
NAME_IMP_3D="implicit-3d"

for name in $NAME_EXP $NAME_IMP_Z $NAME_IMP_3D; do
  echo "INFO: Running CaNS, case $name"
  mv input-${name}.nml input.nml
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
  mkdir -p data-$name
  mv data/* data-$name
done

cp $TESTDIR/*.* .
echo "INFO: Comparing the two runs"
pytest test.py
rm -rf $RUNDIR/data/*.* data-*
