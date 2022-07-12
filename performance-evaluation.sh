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



build_parflow() {
    BP_CMD="make -j16"

    if $BP_CMD ; then
        echo "Build successful."
    else
        echo "Build failed."
        exit 2
    fi
}

# Function Definitions
build_chapel_modules() {
    BCM_FAST=0
    BCM_CALL_ONLY=0
    while test $# -gt 0
    do
        case "$1" in
            --fast) BCM_FAST=1;
                ;;
            --call-only) BCM_CALL_ONLY=1;
                ;;
        esac
        shift
    done

    cd $PARFLOW_DIR/pfsimulator/chapel

    BCM_INCLUDES="-I../parflow_lib -I../parflow_lib \
    -I ../amps/mpi1/ -I../../build/include/ -I../amps/common"

    BCM_COMMAND="chpl $BCM_INCLUDES"
    if [[ $BCM_FAST -eq 1 ]]; then
        BCM_COMMAND="$BCM_COMMAND --fast "
    fi

    if [[ $BCM_CALL_ONLY -eq 1 ]]; then
        BCM_COMMAND="$BCM_COMMAND --scall_only=1 "
    fi

    BCM_COMMAND="$BCM_COMMAND --library chapel_impl.chpl "

    if $BCM_COMMAND; then
        echo "Chapel modules built successfully."
    else
        echo "Chapel module failed to build."
        exit
    fi
}

#Runs the example and saves the execution times. Name of the variant is 
# the first and only argument to this function
run_example() {
RE_VARIANT_NAME=$1

cd $PARFLOW_DIR/test/tcl/washita/tcl_scripts

echo -n "${RE_VARIANT_NAME}," >> $RESULTS_FILE

for i in $RUNSEQ
do
    rm Outputs/cp.out
    echo "${RE_VARIANT_NAME}, run $i"
    tclsh $LW

    RESULT_LINE="$(grep "Total Run Time" Outputs/LW.out.log)"
    RESULT_ARRAY=($RESULT_LINE)
    EXEC_TIME=${RESULT_ARRAY[3]}
    echo $EXEC_TIME
    echo -n " $EXEC_TIME," >> $RESULTS_FILE
done
echo >> $RESULTS_FILE


}

check_correctness() {
    cd $PARFLOW_DIR/test/tcl/washita/tcl_scripts
    diff Outputs/cp.out correct_output/cp.out >> correctness_diff 
    if [ -s correctness_diff ]; then
        echo "Correctness check failed."
    fi
}


LAUNCH_DIR=$(pwd)

CMAKE_DEFS="-DPARFLOW_HAVE_CLM=ON -DCMAKE_INSTALL_PREFIX=${PARFLOW_DIR} \
-DSILO_ROOT=$HOME/silo -DHYPRE_ROOT=$HOME/hypre/src/hypre \
-DPARFLOW_AMPS_LAYER=mpi1 -DCMAKE_BUILD_TYPE=RELWITHDEBINFO \
-DPARFLOW_ENABLE_TIMING=true"

CMAKE_ORIGINAL="cmake ${CMAKE_DEFS} .."

CMAKE_CHAPEL="${CMAKE_ORIGINAL} -DPARFLOW_HAVE_CHAPEL=ON"


echo "Reading command line arguments..."

LW="LW_Tiny.tcl"
RUN_ORIGINAL=0
RUN_CHAPEL=0
RUN_CALL_ONLY=0
RUN_CHAPEL_FAST=0
APEEND=0
RESULT_NAME=performance-evaluation.results
NUMRUNS=2
while test $# -gt 0
do
    case "$1" in
        --original) echo "Running original variant";
            RUN_ORIGINAL=1
            ;;
        --chapel) echo "Running chapel variant";
            RUN_CHAPEL=1
            ;;
        --call-only) echo "Running call only variant";
            RUN_CALL_ONLY=1
            ;;
        --chapel-fast) echo "Running chapel fast variant";
            RUN_CHAPEL_FAST=1
            ;;
        --append) echo "Appending to results file";
            APPEND=1
            ;;
        --outfile) echo "Writing to $2";
            RESULT_NAME=$2
            shift
            ;;
        --n) echo "running $2 times";
            NUMRUNS=$2; shift
            ;;
        *) echo "argument $1"
            ;;
    esac
    shift
done
RUNSEQ=$(seq 1 $NUMRUNS)

echo "Creating results file..."
RESULTS_FILE=$LAUNCH_DIR/$RESULT_NAME
if [[ $APPEND -eq 0 ]]; then
    rm $RESULTS_FILE
fi
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
    run_example original

    
    echo "Done running original variant."
fi


if [[ $RUN_CHAPEL = 1 ]]; then

    echo "Preparing chapel modules..."
    build_chapel_modules


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
    run_example "chapel"
    check_correctness

fi 

if [[ $RUN_CHAPEL_FAST = 1 ]]; then
    echo "Running Chapel --fast Variant"

    echo "Preparing chapel modules..."
    build_chapel_modules --fast


    echo "Clearing build for chapel variant..."
    cd $PARFLOW_DIR/build
    rm -rf ./*

    echo "Configuring..."
    $CMAKE_CHAPEL

    echo "Building..."
    build_parflow

    echo "Installing..."
    make install

    echo "Preparing to run Little Washita Example..."
    run_example "chapel-fast"
    check_correctness
fi

if [[ $RUN_CALL_ONLY -eq 1 ]]; then
    echo "Running Call Only Variant"

    echo "Preparing chapel modules..."
    build_chapel_modules --call-only

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
    run_example "call-only"
fi