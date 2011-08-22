include (FindPackageHandleStandardArgs)

find_path (statgen_INCLUDE_DIRS
  "SamFile.h"
  PATHS
    ${CMAKE_SOURCE_DIR}/../ext/src/statgen/lib/include
)

find_library (statgen_LIBRARIES
  "StatGen"
  PATHS
    ${CMAKE_SOURCE_DIR}/../build/libs/
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(statgen
  "statgen library cannot be found"
  statgen_LIBRARIES
  statgen_INCLUDE_DIRS
)

MARK_AS_ADVANCED(statgen_INCLUDE_DIRS statgen_LIBRARIES)
