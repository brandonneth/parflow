#!/bin/bash

export SILO_DIR=/g/g15/neth3/SparseParFlow/silo-4.11
export HYPRE_DIR=/g/g15/neth3/SparseParFlow/hypre/src/hypre
export PARFLOW_DIR=/g/g15/neth3/SparseParFlow/parflow
configure() {
pushd $PARFLOW_DIR
module load clang/14.0.5
export CXX=clang++
export CC=clang
CMAKE_DEFS="-DCMAKE_C_FLAGS=-g -DPARFLOW_HAVE_CLM=ON -DCMAKE_INSTALL_PREFIX=${PARFLOW_DIR} \
    -DSILO_ROOT=${SILO_DIR} -DHYPRE_ROOT=${HYPRE_DIR} \
    -DPARFLOW_AMPS_LAYER=mpi1 -DCMAKE_BUILD_TYPE=RELWITHDEBINFO  -DPARFLOW_HAVE_OMP=ON"
mkdir build
cd build
cmake $CMAKE_DEFS ..
cd ..
}

build() {
pushd $PARFLOW_DIR/build
make -j
make install
cd ..
 SEARCH="\${MPIEXEC_PREFLAGS} \${PARFLOW_DIR}/bin/\${PROGRAM} \${PFSCRIPT}"
REPLACE="\${MPIEXEC_PREFLAGS} hpcrun -o measurements \${PARFLOW_DIR}/bin/\${PROGRAM} \${PFSCRIPT}"
sed -i 's/PREFLAGS} \${PARFLOW_DIR/PREFLAGS} hpcrun -o toolkit-measurements \${PARFLOW_DIR/g' bin/run
}

run() {
  module load hpctoolkit
  pushd $PARFLOW_DIR/test/tcl/washita/tcl_scripts
  rm -r Outputs/toolkit-measurements
  tclsh LW_Test.tcl
  hpcstruct Outputs/toolkit-measurements
  hpcprof-mpi Outputs/toolkit-measurements -I$PARFLOW_DIR/pfsimulator/+
  
}

