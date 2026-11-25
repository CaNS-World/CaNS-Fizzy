#!/bin/bash
#
# options of the two-fluid method that no other test selects: each is run for a short time on the rising
# bubble of the test `rising_bubble`, and must give a bubble close to that of its reference case, without
# being the same run
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
cp $TESTDIR/../rising_bubble/spheres.in $RUNDIR
run_case() {
  #
  # runs the case named $1, with the lines $2 (if any) added to the `two_fluid` namelist, and the line $3
  # (if any) as the `numerics` namelist
  #
  echo "INFO: Running CaNS, case $1"
  sed -e "s/^nstep.*/nstep = 1000, time_max = 100., tw_max = 0.1/" -e "s/^stop_type.*/stop_type(1:3) = T, F, F/" \
      -e "s/^icheck.*/icheck = 100, iout0d = 100, iout1d = 0, iout2d = 0, iout3d = 0, isave = 1000/" \
      $TESTDIR/../rising_bubble/input.nml > input.nml
  if [ -n "$2" ]; then sed -i "/^&two_fluid/a $2" input.nml; fi
  if [ -n "$3" ]; then printf "\n&numerics\n$3\n/\n" >> input.nml; fi
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
  mkdir -p data-$1
  mv data/* data-$1
}

run_case acdi           ""
run_case acdi-nosdf     "is_sdf_normals = F"
run_case acdi-nobrack   "is_brackbill_normalization = F"
run_case acdi-harmonic  "" "is_harmonic_interpolation = T"
run_case thinc          "ifc_reconstruction = 'thinc'"
run_case thinc-nosdf    "ifc_reconstruction = 'thinc', is_sdf_normals = F"
run_case plic           "ifc_reconstruction = 'plic', ifc_transport = 'sweeps'"
run_case plic-hfnormals "ifc_reconstruction = 'plic', ifc_transport = 'sweeps', is_hf_normals = T"
run_case plic-nolimiter "ifc_reconstruction = 'plic', ifc_transport = 'sweeps'" "is_visc_limiter = F"

cp $TESTDIR/*.* .
echo "INFO: Checking the results"
pytest test.py
rm -rf $RUNDIR/data/*.* data-*
