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

echo "Creating results file..."
RESULTS_FILE=$LAUNCH_DIR/performance-evaluation.results
rm $RESULTS_FILE
touch $RESULTS_FILE

echo "Preparing chapel modules..."

cd $PARFLOW_DIR/pfsimulator/chapel
chpl --library chapel_impl.chpl -I../parflow_lib -I../parflow_lib -I ../amps/mpi1/ -I../../build/include/ -I../amps/common



echo "Original Variant"

cd $PARFLOW_DIR
mkdir build
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
    tclsh LW_Timing.tcl

    RESULT_LINE="$(grep "Total Run Time" Outputs/LW.out.log)"
    RESULT_ARRAY=($RESULT_LINE)
    EXEC_TIME=${RESULT_ARRAY[3]}
    echo $EXEC_TIME
    echo -n " $EXEC_TIME," >> $RESULTS_FILE
done

echo >> $RESULTS_FILE


echo "Done running original variant."


cd $PARFLOW_DIR/build

echo "Clearing build for chapel variant..."
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
    tclsh LW_Timing.tcl

    RESULT_LINE="$(grep "Total Run Time" Outputs/LW.out.log)"
    RESULT_ARRAY=($RESULT_LINE)
    EXEC_TIME=${RESULT_ARRAY[3]}
    echo $EXEC_TIME
    echo -n " $EXEC_TIME," >> $RESULTS_FILE
done
