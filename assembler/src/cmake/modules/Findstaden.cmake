include (FindPackageHandleStandardArgs)

find_path (staden_INCLUDE_DIRS
  "io_lib"
  PATHS
    /usr/local/include
    /usr/include
		${CMAKE_SOURCE_DIR}/../build/ext/staden/include
)

find_library (staden_LIBRARIES
  "staden-read"
  PATHS
    /usr/local/lib
    /usr/lib
		${CMAKE_SOURCE_DIR}/../build/ext/staden/lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(staden
  "staden library cannot be found"
  staden_LIBRARIES
  staden_INCLUDE_DIRS
)

MARK_AS_ADVANCED(staden_INCLUDE_DIRS staden_LIBRARIES)
