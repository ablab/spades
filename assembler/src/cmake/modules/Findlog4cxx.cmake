include (FindPackageHandleStandardArgs)

find_path (log4cxx_INCLUDE_DIRS
  "log4cxx"
  PATHS
    /usr/local/include
    /usr/include
)

find_library (log4cxx_LIBRARIES
  "log4cxx"
  PATHS
    /usr/local/lib
    /usr/lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(log4cxx
  "log4cxx library cannot be found"
  log4cxx_LIBRARIES
  log4cxx_INCLUDE_DIRS
)

MARK_AS_ADVANCED(log4cxx_INCLUDE_DIRS log4cxx_LIBRARIES)
