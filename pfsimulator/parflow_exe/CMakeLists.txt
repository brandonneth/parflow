

add_executable(parflow main.c)

target_include_directories(parflow PRIVATE "../chapel/lib")
target_include_directories(parflow PRIVATE "~/chapel/lib/darwin/clang/x86_64/cpu-native/loc-flat/comm-none/tasks-fifo/tmr-generic/unwind-none/mem-cstdlib/atomics-cstdlib/hwloc-none/re2-none/fs-none/lib_pic-none/san-none")
target_include_directories(parflow PRIVATE "~/chapel/runtime/include")
set (CHPL_HOME /users/brandonneth/chapel)
set (CHPL_RUNTIME_LIB ${CHPL_HOME}/lib)
set (CHPL_RUNTIME_INCL ${CHPL_HOME}/runtime/include)
set (CHPL_THIRD_PARTY ${CHPL_HOME}/third-party)

target_include_directories(parflow PRIVATE ${CHPL_HOME}/modules/internal ${CHPL_RUNTIME_INCL}/localeModels/flat ${CHPL_RUNTIME_INCL}/localeModels ${CHPL_RUNTIME_INCL}/comm/none ${CHPL_RUNTIME_INCL}/comm ${CHPL_RUNTIME_INCL}/tasks/fifo ${CHPL_RUNTIME_INCL}/. ${CHPL_RUNTIME_INCL}/./qio ${CHPL_RUNTIME_INCL}/./atomics/cstdlib ${CHPL_RUNTIME_INCL}/./mem/cstdlib ${CHPL_THIRD_PARTY}/utf8-decoder)
target_link_directories(parflow PRIVATE ${CHPL_RUNTIME_LIB}/darwin/clang/x86_64/cpu-native/loc-flat/comm-none/tasks-fifo/tmr-generic/unwind-none/mem-cstdlib/atomics-cstdlib/hwloc-none/re2-none/fs-none/lib_pic-none/san-none ../chapel/lib)

target_link_libraries(parflow chpl)


target_link_libraries(parflow pfsimulator)

if( ${PARFLOW_HAVE_ETRACE} )
  target_include_directories(parflow PUBLIC "../third_party/etrace")
  target_link_libraries(parflow etrace)
endif( ${PARFLOW_HAVE_ETRACE} )

if( ${PARFLOW_HAVE_RMM} )
  target_link_libraries(parflow rmm)
  target_include_directories(parflow PRIVATE ${RMM_INCLUDE})
  target_compile_definitions(parflow PRIVATE PARFLOW_HAVE_RMM)
endif( ${PARFLOW_HAVE_RMM} )

if( ${PARFLOW_HAVE_CLM} )
  target_link_libraries(parflow pfclm)
endif( ${PARFLOW_HAVE_CLM} )

if (${PARFLOW_HAVE_HYPRE})
  target_link_libraries (parflow ${HYPRE_LIBRARIES})
endif (${PARFLOW_HAVE_HYPRE})

if (${PARFLOW_HAVE_MPI})
  # In CMake 3.13 this could be target_link_options
  target_link_libraries(parflow ${MPI_LINK_FLAGS})
  target_link_libraries (parflow ${MPI_LIBRARIES})
endif (${PARFLOW_HAVE_MPI})

if (${PARFLOW_HAVE_SILO})
  target_link_libraries (parflow ${SILO_LIBRARIES})
endif (${PARFLOW_HAVE_SILO})

if (${PARFLOW_HAVE_NETCDF})
  target_link_libraries (parflow ${NETCDF_LIBRARIES})

  target_link_libraries (parflow ${CURL_LIBRARIES})
endif (${PARFLOW_HAVE_NETCDF})

if (${PARFLOW_HAVE_HDF5})
  target_link_libraries (parflow ${HDF5_LIBRARIES})

  if (${PARFLOW_HAVE_NETCDF})
    target_link_libraries (parflow ${HDF5_HL_LIBRARIES})

  endif (${PARFLOW_HAVE_NETCDF})
endif (${PARFLOW_HAVE_HDF5})

if (${PARFLOW_HAVE_ZLIB})
  target_link_libraries (parflow ${ZLIB_LIBRARIES})
endif (${PARFLOW_HAVE_ZLIB})

if (${PARFLOW_HAVE_SZLIB})
  target_link_libraries (parflow ${SZLIB_LIBRARIES})
endif (${PARFLOW_HAVE_SZLIB})

if (${PARFLOW_HAVE_SLURM})
  target_link_libraries (parflow ${SLURM_LIBRARIES})
endif (${PARFLOW_HAVE_SLURM})

if( ${PARFLOW_ENABLE_PROFILING} )
  set_target_properties(parflow PROPERTIES LINK_FLAGS ${PARFLOW_PROFILE_OPTS})
endif( ${PARFLOW_ENABLE_PROFILING} )

if ( DEFINED PARFLOW_LINKER_FLAGS)
   set_target_properties(parflow PROPERTIES LINK_FLAGS ${PARFLOW_LINKER_FLAGS})
endif ( DEFINED PARFLOW_LINKER_FLAGS)


install(TARGETS parflow DESTINATION bin)
