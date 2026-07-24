#!/bin/bash
set -eu

TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
RUNDIR=$CANSDIR/run

rm -rf "$RUNDIR"
cp "$TESTDIR/build.conf" "$CANSDIR/build.conf"
cd "$CANSDIR"
make allclean
make libs
make -j

CHANNELDIR=$RUNDIR/zero_net_channel
mkdir -p "$CHANNELDIR/data"
cp "$RUNDIR/cans" "$CHANNELDIR"
cp "$TESTDIR/input.nml" "$CHANNELDIR"
cd "$CHANNELDIR"
mpirun -n 4 --oversubscribe ./cans
cd data
pytest "$TESTDIR/test.py"

CAPILLARYDIR=$RUNDIR/zero_net_capillary
mkdir -p "$CAPILLARYDIR/data"
cp "$RUNDIR/cans" "$CAPILLARYDIR"
cp "$TESTDIR/input_capillary.nml" "$CAPILLARYDIR/input.nml"
cp "$TESTDIR/spheres.in" "$CAPILLARYDIR"
cd "$CAPILLARYDIR"
mpirun -n 4 --oversubscribe ./cans
cd data
pytest "$TESTDIR/test_capillary.py"

CAPILLARYEQDIR=$RUNDIR/zero_net_capillary_equal_density
mkdir -p "$CAPILLARYEQDIR/data"
cp "$RUNDIR/cans" "$CAPILLARYEQDIR"
cp "$TESTDIR/input_capillary.nml" "$CAPILLARYEQDIR/input.nml"
cp "$TESTDIR/spheres.in" "$CAPILLARYEQDIR"
sed -i 's/rho12(1:2) = 10., 1./rho12(1:2) = 1., 1./' "$CAPILLARYEQDIR/input.nml"
cd "$CAPILLARYEQDIR"
mpirun -n 4 --oversubscribe ./cans
cd data
pytest "$TESTDIR/test_capillary_momentum.py"

BOUSSINESQDIR=$RUNDIR/zero_net_boussinesq
mkdir -p "$BOUSSINESQDIR/data"
cp "$RUNDIR/cans" "$BOUSSINESQDIR"
cp "$TESTDIR/input_boussinesq.nml" "$BOUSSINESQDIR/input.nml"
cd "$BOUSSINESQDIR"
mpirun -n 4 --oversubscribe ./cans
cd data
pytest "$TESTDIR/test_boussinesq.py"

INVALIDDIR=$RUNDIR/zero_net_invalid
mkdir -p "$INVALIDDIR"
for CASE in conflict nonperiodic open; do
  mkdir -p "$INVALIDDIR/$CASE/data"
  cp "$RUNDIR/cans" "$INVALIDDIR/$CASE"
  cp "$TESTDIR/input.nml" "$INVALIDDIR/$CASE/input.nml"
  sed -i 's/dims(1:2) = 2, 2/dims(1:2) = 1, 1/' "$INVALIDDIR/$CASE/input.nml"
done

sed -i 's/bforce(1:3) = 0., 0., 0./bforce(1:3) = 1., 0., 0./' "$INVALIDDIR/conflict/input.nml"
if (cd "$INVALIDDIR/conflict" && mpirun -n 1 ./cans >run.log 2>&1); then
  echo "Conflicting fixed and dynamic forcing was not rejected."
  exit 1
fi
grep -q 'bforce must be zero' "$INVALIDDIR/conflict/run.log"

sed -i "s/'P','P',  'P','P',  'D','D'/'D','D',  'P','P',  'D','D'/" "$INVALIDDIR/nonperiodic/input.nml"
sed -i "s/'P','P',  'P','P',  'N','N'/'N','N',  'P','P',  'N','N'/" "$INVALIDDIR/nonperiodic/input.nml"
if (cd "$INVALIDDIR/nonperiodic" && mpirun -n 1 ./cans >run.log 2>&1); then
  echo "Forcing in a nonperiodic direction was not rejected."
  exit 1
fi
grep -q 'requires periodic pressure BCs' "$INVALIDDIR/nonperiodic/run.log"

sed -i \
  's/bcvel(0:1,1:3,3) = 0.,0.,  0.,0.,  0.,0./bcvel(0:1,1:3,3) = 0.,0.,  0.,0.,  1.,1./' \
  "$INVALIDDIR/open/input.nml"
if (cd "$INVALIDDIR/open" && mpirun -n 1 ./cans >run.log 2>&1); then
  echo "An open transverse boundary was not rejected."
  exit 1
fi
grep -q 'periodic or impermeable transverse boundaries' "$INVALIDDIR/open/run.log"
