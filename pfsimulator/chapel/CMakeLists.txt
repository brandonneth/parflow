#Definitions for linking chapel runtime and this module with parflow
set (CHPL_HOME $ENV{CHPL_HOME})
set (CHPL_RUNTIME_LIB ${CHPL_HOME}/lib)
set (CHPL_RUNTIME_INCL ${CHPL_HOME}/runtime/include)
set (CHPL_THIRD_PARTY ${CHPL_HOME}/third-party)

set (CHPL_INCLUDE_DIRECTORIES ${CHPL_HOME}/modules/internal ${CHPL_RUNTIME_INCL}/localeModels/flat $(CHPL_RUNTIME_INCL)/localeModels ${CHPL_RUNTIME_INCL}/comm/none ${CHPL_RUNTIME_INCL}/comm ${CHPL_RUNTIME_INCL}/tasks/qthreads ${CHPL_RUNTIME_INCL}/. ${CHPL_RUNTIME_INCL}/./qio ${CHPL_RUNTIME_INCL}/./atomics/cstdlib ${CHPL_RUNTIME_INCL}/./mem/jemalloc ${CHPL_THIRD_PARTY}/utf8-decoder -DCHPL_JEMALLOC_PREFIX=chpl_je_ ${CHPL_THIRD_PARTY}/gmp/install/cray-xc-x86_64-none-cray-prgenv-cray-none/include ${CHPL_THIRD_PARTY}/hwloc/install/cray-xc-x86_64-none-cray-prgenv-cray-none-flat/include ${CHPL_THIRD_PARTY}/qthread/install/cray-xc-x86_64-none-cray-prgenv-cray-none-flat-jemalloc-bundled/include ${CHPL_THIRD_PARTY}/jemalloc/install/target/cray-xc-x86_64-none-cray-prgenv-cray-none/include ${CHPL_THIRD_PARTY}/re2/install/cray-xc-x86_64-none-cray-prgenv-cray-none/include ${CMAKE_CURRENT_LIST_DIR}/lib)
set (CHPL_LINK_DIRECTORIES ${CHPL_RUNTIME_LIB}/cray-xc/cray-prgenv-cray/x86_64/cpu-none/loc-flat/comm-none/tasks-qthreads/tmr-generic/unwind-none/mem-jemalloc/atomics-cstdlib/hwloc-bundled/re2-bundled/fs-none/lib_pic-none/san-none  ${CHPL_THIRD_PARTY}/gmp/install/cray-xc-x86_64-none-cray-prgenv-cray-none/lib $(CHPL_THIRD_PARTY)/gmp/install/cray-xc-x86_64-none-cray-prgenv-cray-none/lib ${CHPL_THIRD_PARTY}/hwloc/install/cray-xc-x86_64-none-cray-prgenv-cray-none-flat/lib ${CHPL_THIRD_PARTY}/hwloc/install/cray-xc-x86_64-none-cray-prgenv-cray-none-flat/lib ${CHPL_THIRD_PARTY}/qthread/install/cray-xc-x86_64-none-cray-prgenv-cray-none-flat-jemalloc-bundled/lib ${CHPL_THIRD_PARTY}/qthread/install/cray-xc-x86_64-none-cray-prgenv-cray-none-flat-jemalloc-bundled/lib ${CHPL_THIRD_PARTY}/hwloc/install/cray-xc-x86_64-none-cray-prgenv-cray-none-flat/lib  ${CHPL_THIRD_PARTY}/jemalloc/install/target/cray-xc-x86_64-none-cray-prgenv-cray-none/lib  ${CHPL_THIRD_PARTY}/re2/install/cray-xc-x86_64-none-cray-prgenv-cray-none/lib  ${CHPL_THIRD_PARTY}/re2/install/cray-xc-x86_64-none-cray-prgenv-cray-none/lib ${CMAKE_CURRENT_LIST_DIR}/lib)

set (CHPL_LINK_LIBRARIES chpl gmp hwloc jemalloc re2 m numa pthread qthread chapel_impl)

#for a target named tgt, add the following lines:
#include(<Path to this file>)
#target_include_directories(tgt PRIVATE $CHPL_INCLUDE_DIRECTORIES)
#target_link_directories(tgt PRIVATE $CHPL_LINK_DIRECTORIES)
#target_link_libraries(tgt $CHPL_LINK_LIBRARIES)
