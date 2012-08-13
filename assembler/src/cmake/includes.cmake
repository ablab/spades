# -*- cmake -*-

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${SPADES_MAIN_INCLUDE_DIR} ${SPADES_BUILT_INCLUDE_DIR})
include_directories("${EXT_DIR}/include")
include_directories("${ZLIB_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIRS}")

if (SPADES_USE_TCMALLOC)
  include_directories("${GOOGLE_PERFTOOLS_INCLUDE_DIR}")
endif()
