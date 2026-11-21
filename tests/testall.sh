#!/bin/bash
#
# builds the variants of CaNS-Fizzy that the tests need, and runs the tests of each variant;
# the variants and their tests are the same as in the CI (.github/workflows/run-compiler.yml)
#
# usage, from this directory: ./testall.sh [variant ...]     (default: all variants)
#
# optional environment variables:
#   BUILD_FLAGS: build flags added to those of the variant, e.g., BUILD_FLAGS="FFLAGS_DEBUG_MAX=1 FFLAGS_OPT=0"
#   TESTS      : runs only the tests of this list, e.g., TESTS="rising_bubble restart_two_fluid"
#   NJOBS      : number of parallel compilation jobs (default: 4)
#
# the dependencies are built if needed (`make libs`); the compiler is the one of `build.conf`,
# and the `hypre` variant needs the HYPRE library (see `docs/INFO_COMPILING.md`);
# the build of the last variant is left in place; the output of each test is saved in
# `run/testall/`, and the exit status is non-zero if a build or a test fails
#
TESTSDIR=$(cd "$(dirname "$0")" && pwd)
CANSDIR=$TESTSDIR/..
LOGDIR=$CANSDIR/run/testall
#
# build flags and tests of each variant
#
VARIANTS="acdi thinc balanced hypre plic"
declare -A FLAGS TESTS_OF
FLAGS[acdi]="CONSTANT_COEFFS_POISSON=1 INTERFACE_CAPTURING_VOF=0 BALANCED_CAPILLARY_PRESSURE_SPLIT=0"
FLAGS[thinc]="CONSTANT_COEFFS_POISSON=1 INTERFACE_CAPTURING_VOF=1 BALANCED_CAPILLARY_PRESSURE_SPLIT=0"
FLAGS[balanced]="CONSTANT_COEFFS_POISSON=1 INTERFACE_CAPTURING_VOF=0 BALANCED_CAPILLARY_PRESSURE_SPLIT=1"
FLAGS[hypre]="CONSTANT_COEFFS_POISSON=0 INTERFACE_CAPTURING_VOF=0"
FLAGS[plic]="CONSTANT_COEFFS_POISSON=1 INTERFACE_CAPTURING_SPLIT_VOF=1"
TESTS_OF[acdi]="lid_driven_cavity differentially_heated_cavity deforming_circle rising_bubble \
                start_stop restart_two_fluid implicit_diffusion poisson taylor_green_vortex_2d"
TESTS_OF[thinc]="deforming_circle rising_bubble restart_two_fluid implicit_diffusion"
TESTS_OF[balanced]="zero_net_acceleration rising_bubble restart_two_fluid"
TESTS_OF[hypre]="rising_bubble restart_two_fluid implicit_diffusion"
TESTS_OF[plic]="plic_unit plic_advection plic_momentum viscosity_limiter rising_bubble restart_two_fluid"
#
if [ $# -gt 0 ]; then VARIANTS="$*"; fi
for variant in $VARIANTS; do
  if [ -z "${FLAGS[$variant]+x}" ]; then echo "ERROR: unknown variant '$variant' (available: ${!FLAGS[*]})"; exit 1; fi
done
#
mkdir -p $LOGDIR
cd $CANSDIR
echo "INFO: Building the dependencies, if needed"
make libs > $LOGDIR/build-libs.log 2>&1 || { echo "ERROR: 'make libs' failed, see $LOGDIR/build-libs.log"; exit 1; }
nfail=0
summary=""
for variant in $VARIANTS; do
  echo "INFO: Building variant $variant"
  make clean > /dev/null 2>&1
  if ! make -j ${NJOBS:-4} ${FLAGS[$variant]} $BUILD_FLAGS > $LOGDIR/build-$variant.log 2>&1; then
    echo "ERROR: the build of variant $variant failed, see $LOGDIR/build-$variant.log"
    summary+=$(printf '%-10s %-30s %s' $variant '(build)' FAILED)$'\n'
    nfail=$((nfail+1))
    continue
  fi
  for test in ${TESTS_OF[$variant]}; do
    if [ -n "$TESTS" ] && [[ " $TESTS " != *" $test "* ]]; then continue; fi
    echo "INFO: Running test $test ($variant)"
    tbeg=$SECONDS
    if (cd $TESTSDIR/$test && bash ./testit.sh) > $LOGDIR/$variant-$test.log 2>&1; then result=passed; else result=FAILED; nfail=$((nfail+1)); fi
    summary+=$(printf '%-10s %-30s %s %5d s' $variant $test $result $((SECONDS-tbeg)))$'\n'
  done
done
echo
echo "Summary (output of each test in $LOGDIR):"
echo -n "$summary"
if [ $nfail -ne 0 ]; then
  echo "ERROR: $nfail failure(s)"
  exit 1
fi
echo "INFO: all tests passed"
