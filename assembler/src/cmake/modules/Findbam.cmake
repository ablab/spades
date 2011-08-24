include (FindPackageHandleStandardArgs)

find_path (bam_INCLUDE_DIRS
  "bam.h"
  PATHS
    ${CMAKE_SOURCE_DIR}/../ext/include/statgen/samtools
)

find_library (bam_LIBRARIES
  "bam"
  PATHS
    ${CMAKE_SOURCE_DIR}/../build/ext/lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(bam
  "bam library cannot be found"
  bam_LIBRARIES
  bam_INCLUDE_DIRS
)

MARK_AS_ADVANCED(bam_INCLUDE_DIRS bam_LIBRARIES)
