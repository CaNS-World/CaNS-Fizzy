#!/bin/bash
#
# targets of the forcing: the zero-net-heating forcing of a scalar (`is_sforced`, `scalf`), and the bulk
# velocity of the initial condition along a forced direction (`velf`)
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
cp $TESTDIR/*.* .
#
# conduction : a single fluid between two cold walls, held at its mean by the forcing, against the exact steady state
# bubble     : a rising bubble with other properties than the fluid around it: the total of rho*cp*s is kept
# channel    : a bubble in a forced channel flow: the mass-averaged velocity is the target, and is kept
#
# the first two also with implicit diffusion: in all directions at five times the time step (`-impdiff`), and
# along z (`-impdiff-z`)
#
for name in conduction conduction-impdiff bubble bubble-impdiff-z channel; do
  base=${name%%-*}
  cp $TESTDIR/input-$base.nml input.nml
  if [ -f $TESTDIR/spheres-$base.in ]; then cp $TESTDIR/spheres-$base.in spheres.in; fi
  case $name in
    conduction-impdiff) sed -i "s/dt_f = 4.e-3/dt_f = 2.e-2/" input.nml; printf "\n&numerics\nis_impdiff = T\n/\n" >> input.nml;;
    bubble-impdiff-z)   sed -i "s/^dims(1:2) *=.*/dims(1:2) = 1, 4, ipencil_axis = 3/" input.nml
                        printf "\n&numerics\nis_impdiff = T, is_impdiff_1d = T\n/\n" >> input.nml;;
  esac
  #
  # interface method of the variant under test, if set (see `tests/testall.sh`): a line for the `two_fluid` namelist
  #
  if [ -n "$TEST_TWO_FLUID_INPUT" ] && [ $base != conduction ]; then sed -i "/^&two_fluid/a $TEST_TWO_FLUID_INPUT" input.nml; fi
  echo "INFO: Running CaNS, case $name"
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
  echo "INFO: Checking the results, case $name"
  FORCING_CASE=$base pytest test.py
  rm -rf $RUNDIR/data/*.* spheres.in
done
