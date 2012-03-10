include (FindPackageHandleStandardArgs)

find_path (python_INCLUDE_DIRS
  "Python.h"
  PATHS
    /usr/local/include/python2.6
    /usr/include/python2.6
)

find_library (python_LIBRARIES
  "python2.6"
  PATHS
    /usr/local/lib
    /usr/lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(python
  "python library cannot be found"
  python_LIBRARIES
  python_INCLUDE_DIRS
)

MARK_AS_ADVANCED(python_INCLUDE_DIRS python_LIBRARIES)
