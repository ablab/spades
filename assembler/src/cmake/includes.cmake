# -*- cmake -*-

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem ")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
include_directories(${SPADES_MAIN_INCLUDE_DIR} ${SPADES_BUILT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR} ${SPADES_MODULES_DIR}) 
include_directories(SYSTEM "${EXT_DIR}/include")
include_directories(SYSTEM "${ZLIB_INCLUDE_DIRS}")
include_directories(SYSTEM "${Boost_INCLUDE_DIRS}")

if (SPADES_USE_TCMALLOC)
  include_directories("${GOOGLE_PERFTOOLS_INCLUDE_DIR}")
endif()
