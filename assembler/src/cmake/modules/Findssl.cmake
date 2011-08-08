include (FindPackageHandleStandardArgs)

find_path (ssl_INCLUDE_DIRS
  "openssl"
  PATHS
    /usr/local/include
    /usr/include
)

find_library (ssl_LIBRARIES
  "ssl"
  PATHS
    /usr/local/lib
    /usr/lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(ssl
  "ssl library cannot be found"
  ssl_LIBRARIES
  ssl_INCLUDE_DIRS
)

MARK_AS_ADVANCED(ssl_INCLUDE_DIRS ssl_LIBRARIES)
