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

for name in one same split; do
  echo "INFO: Running CaNS, case $name"
  cp $TESTDIR/input.nml input.nml
  cp spheres-same.in spheres.in
  if [ $name = one   ]; then sed -i "s/^nmarker = 2/nmarker = 1/" input.nml; fi
  if [ $name = split ]; then cp spheres-split.in spheres.in; fi
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
echo "INFO: Checking the results"
pytest test.py
rm -rf $RUNDIR/data/*.* data-*
