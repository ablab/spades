include (FindPackageHandleStandardArgs)

find_path (fftw3_INCLUDE_DIRS
  "fftw3.h"
  PATHS
    /usr/local/include
    /usr/include
)

find_library (fftw3_LIBRARIES
  "fftw3"
  PATHS
    /usr/local/lib
    /usr/lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(fftw3
  "fftw3 library cannot be found"
  fftw3_LIBRARIES
  fftw3_INCLUDE_DIRS
)

MARK_AS_ADVANCED(fftw3_INCLUDE_DIRS fftw3_LIBRARIES)
