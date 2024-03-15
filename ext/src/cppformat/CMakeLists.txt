project(FORMAT)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include)

set(FMT_SOURCES format.cc)

include(CheckSymbolExists)
check_symbol_exists(open fcntl.h HAVE_OPEN)
if (HAVE_OPEN)
  add_definitions(-DFMT_USE_FILE_DESCRIPTORS=1)
  set(FMT_SOURCES ${FMT_SOURCES} posix.cc)
endif ()

add_library(format ${FMT_SOURCES})
set_target_properties(format PROPERTIES COMPILE_FLAGS  "-Wall -Wextra -pedantic")
