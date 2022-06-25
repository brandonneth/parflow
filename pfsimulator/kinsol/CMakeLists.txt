
set (SRC_FILES iterativ.c kinsol.c kinspgmr.c llnlmath.c spgmr.c)

add_library(pfkinsol ${SRC_FILES})

if(${PARFLOW_HAVE_CUDA})
  target_include_directories(pfkinsol PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif(${PARFLOW_HAVE_CUDA})

if( ${PARFLOW_HAVE_RMM} )
  target_link_libraries(pfkinsol rmm)
  target_include_directories(pfkinsol PRIVATE ${RMM_INCLUDE})
  target_compile_definitions(pfkinsol PRIVATE PARFLOW_HAVE_RMM)
endif( ${PARFLOW_HAVE_RMM} )

if (${PARFLOW_HAVE_MPI})
  target_include_directories (pfkinsol PUBLIC "${MPI_C_INCLUDE_PATH}")
endif (${PARFLOW_HAVE_MPI})

set (CHPL_HOME /home/users/neth/chapel)
set (CHPL_RUNTIME_LIB ${CHPL_HOME}/lib)
set (CHPL_RUNTIME_INCL ${CHPL_HOME}/runtime/include)
set (CHPL_THIRD_PARTY ${CHPL_HOME}/third-party)

target_include_directories(pfkinsol PRIVATE /home/users/neth/chapel/modules/internal ${CHPL_RUNTIME_INCL}/localeModels/flat ${CHPL_RUNTIME_INCL}/localeModels ${CHPL_RUNTIME_INCL}/comm/none ${CHPL_RUNTIME_INCL}/comm ${CHPL_RUNTIME_INCL}/tasks/fifo ${CHPL_RUNTIME_INCL}/. ${CHPL_RUNTIME_INCL}/./qio ${CHPL_RUNTIME_INCL}/./atomics/cstdlib ${CHPL_RUNTIME_INCL}/./mem/cstdlib ${CHPL_THIRD_PARTY}/utf8-decoder ../chapel/lib)

target_link_directories(pfkinsol PRIVATE ${CHPL_RUNTIME_LIB}/cray-xc/cray-prgenv-cray/x86_64/cpu-none/loc-flat/comm-none/tasks-fifo/tmr-generic/unwind-none/mem-cstdlib/atomics-cstdlib/hwloc-none/re2-none/fs-none/lib_pic-none/san-none ../chapel/lib)

target_link_libraries(pfkinsol chpl)

install (TARGETS pfkinsol DESTINATION lib)
