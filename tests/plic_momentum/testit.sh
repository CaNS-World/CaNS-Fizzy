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
cp $TESTDIR/input-*.nml $TESTDIR/spheres-*.in $RUNDIR

for name in drop drop-thinc restart-one-step restart-two-step-1 restart-two-step-2; do
  echo "INFO: Running CaNS, case $name"
  mv input-${name}.nml input.nml
  cp spheres-$(echo $name | sed 's/-.*//g').in spheres.in
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
  dirname="data-$(echo $name | sed 's/-[1-2]$//g')/"
  mkdir -p $dirname
  if [[ $name == *"two-step"* ]]; then
    cp -r data/* $dirname
  else
    mv data/* $dirname
  fi
done

cp $TESTDIR/*.* .
echo "INFO: Checking the results"
pytest test.py
rm -rf $RUNDIR/data/*.* data-* spheres-*.in
