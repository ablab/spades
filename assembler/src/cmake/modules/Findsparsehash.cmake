include (FindPackageHandleStandardArgs)

find_path (sparsehash_INCLUDE_DIRS
  "sparse_hash_map"
  PATHS
    /usr/local/include
    /usr/include
    ${CMAKE_SOURCE_DIR}/../ext/include/google
)

#find_library (sparsehash_LIBRARIES
#  "sparsehash"
#  PATHS
#    /usr/local/lib
#    /usr/lib
#    ${CMAKE_SOURCE_DIR}/../build/ext/sparsehash/lib
#)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(sparsehash
  "sparsehash library cannot be found"
#  sparsehash_LIBRARIES
  sparsehash_INCLUDE_DIRS
)

MARK_AS_ADVANCED(sparsehash_INCLUDE_DIRS sparsehash_LIBRARIES)
