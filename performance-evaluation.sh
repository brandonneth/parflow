#!/bin/bash

# Script to compare performance of the original parflow against the chapel
#  interop variant
#
# Arguments: Number of repetitions to run, default 5.
#
# Assumptions:
#   PARFLOW_DIR is set to the parflow base directory,
#   CHPL_HOME is set to the chapel base directory
#   Others to come
#
# Outputs:
#   performance-evaluation.results: execution times for original and chapel
#       variants. Each line contains the variant ("original" or "chapel") and
#       a comma separated list of the execution times. Example: original, 1.2, 1.3

LAUNCH_DIR=$(pwd)

CMAKE_DEFS="-DPARFLOW_HAVE_CLM=ON -DCMAKE_INSTALL_PREFIX=${PARFLOW_DIR} \
-DSILO_ROOT=$HOME/silo -DHYPRE_ROOT=$HOME/hypre/src/hypre \
-DPARFLOW_AMPS_LAYER=mpi1 -DCMAKE_BUILD_TYPE=RELWITHDEBINFO \
-DPARFLOW_ENABLE_TIMING=true"

CMAKE_ORIGINAL="cmake ${CMAKE_DEFS} .."

CMAKE_CHAPEL="${CMAKE_ORIGINAL} -DPARFLOW_HAVE_CHAPEL=ON"
RUNSEQ=$(seq 1 5)

echo "Reading command line arguments..."

LW="LW_Test.tcl"
RUN_ORIGINAL=1
RUN_CHAPEL=1
RUN_CALL_ONLY=1
RUN_CHAPEL_FAST=1
while test $# -gt 0
do
    case "$1" in
        --no-original) echo "Not running original variant";
            RUN_ORIGINAL=0
            ;;
        --no-chapel) echo "Not running chapel variant";
            RUN_CHAPEL=0
            ;;
        --no-call-only) echo "Not running call only variant";
            RUN_CALL_ONLY=0
            ;;
        --no-chapel-fast) echo "Not running chapel fast variant";
            RUN_CHAPEL_FAST=0
            ;;
        *) echo "argument $1"
            ;;
    esac
    shift
done


echo "Creating results file..."
RESULTS_FILE=$LAUNCH_DIR/performance-evaluation.results
touch $RESULTS_FILE

echo "Creating build directory..."
cd $PARFLOW_DIR
mkdir build

if [[ $RUN_ORIGINAL -eq 1 ]]; then
cd build
echo "Clearing build directory..."
rm -rf ./*

echo "Configuring..."
$CMAKE_ORIGINAL

echo "Building..."
make -j8

echo "Installing..."
make install

echo "Preparing to run Little Washita Example..."
cd $PARFLOW_DIR/test/tcl/washita/tcl_scripts



echo -n "original," >> $RESULTS_FILE
for i in $RUNSEQ
do
    echo "Original, run $i"
    tclsh $LW

    RESULT_LINE="$(grep "Total Run Time" Outputs/LW.out.log)"
    RESULT_ARRAY=($RESULT_LINE)
    EXEC_TIME=${RESULT_ARRAY[3]}
    echo $EXEC_TIME
    echo -n " $EXEC_TIME," >> $RESULTS_FILE
done
echo >> $RESULTS_FILE

echo "Done running original variant."
fi


if [[ $RUN_CHAPEL = 1 ]]; then

echo "Preparing chapel modules..."

cd $PARFLOW_DIR/pfsimulator/chapel
chpl --library chapel_impl.chpl -I../parflow_lib -I../parflow_lib \
-I ../amps/mpi1/ -I../../build/include/ -I../amps/common


echo "Clearing build for chapel variant..."
cd $PARFLOW_DIR/build
rm -rf ./*

echo "Configuring..."
$CMAKE_CHAPEL

echo "Building..."
make -j8

echo "Installing..."
make install

echo "Preparing to run Little Washita Example..."
cd $PARFLOW_DIR/test/tcl/washita/tcl_scripts

echo -n "chapel," >> $RESULTS_FILE
for i in $RUNSEQ
do
    echo "Chapel, run $i"
    tclsh $LW

    RESULT_LINE="$(grep "Total Run Time" Outputs/LW.out.log)"
    RESULT_ARRAY=($RESULT_LINE)
    EXEC_TIME=${RESULT_ARRAY[3]}
    echo $EXEC_TIME
    echo -n " $EXEC_TIME," >> $RESULTS_FILE
done
echo >> $RESULTS_FILE
fi

if [[ $RUN_CHAPEL_FAST = 1 ]]; then

echo "Preparing chapel modules..."

cd $PARFLOW_DIR/pfsimulator/chapel
chpl --fast --library chapel_impl.chpl -I../parflow_lib -I../parflow_lib \
-I ../amps/mpi1/ -I../../build/include/ -I../amps/common


echo "Clearing build for chapel variant..."
cd $PARFLOW_DIR/build
rm -rf ./*

echo "Configuring..."
$CMAKE_CHAPEL

echo "Building..."
make -j8

echo "Installing..."
make install

echo "Preparing to run Little Washita Example..."
cd $PARFLOW_DIR/test/tcl/washita/tcl_scripts

echo -n "chapel-fast," >> $RESULTS_FILE
for i in $RUNSEQ
do
    echo "Chapel Fast, run $i"
    tclsh $LW

    RESULT_LINE="$(grep "Total Run Time" Outputs/LW.out.log)"
    RESULT_ARRAY=($RESULT_LINE)
    EXEC_TIME=${RESULT_ARRAY[3]}
    echo $EXEC_TIME
    echo -n " $EXEC_TIME," >> $RESULTS_FILE
done
echo >> $RESULTS_FILE
fi

if [[ $RUN_CALL_ONLY -eq 1 ]]; then
echo "Running Call Only Variant"

echo "Preparing chapel modules..."
cd $PARFLOW_DIR/pfsimulator/chapel
chpl --fast --library chapel_impl.chpl -I../parflow_lib -I../parflow_lib \
-I ../amps/mpi1/ -I../../build/include/ -I../amps/common -scall_only=1

cd $PARFLOW_DIR/build
echo "Clearing build directory..."
rm -rf ./*

echo "Configuring..."
$CMAKE_CHAPEL -DPARFLOW_CALL_ONLY=On

echo "Building..."
make -j8

echo "Installing..."
make install

echo "Preparing to run Little Washita Example..."
cd $PARFLOW_DIR/test/tcl/washita/tcl_scripts

echo -n "call-only," >> $RESULTS_FILE
for i in $RUNSEQ
do
    echo "Chapel call only, run $i"
    tclsh $LW

    RESULT_LINE="$(grep "Total Run Time" Outputs/LW.out.log)"
    RESULT_ARRAY=($RESULT_LINE)
    EXEC_TIME=${RESULT_ARRAY[3]}
    echo $EXEC_TIME
    echo -n " $EXEC_TIME," >> $RESULTS_FILE
done
echo >> $RESULTS_FILE
fi
