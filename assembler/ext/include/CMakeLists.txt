# -*- cmake -*-

set(BLAZE_LIBRARY_PATH ${EXT_DIR}/include/blaze)
list(APPEND CMAKE_MODULE_PATH ${BLAZE_LIBRARY_PATH}/cmake)
include(Blaze_Import)

Blaze_Import(
   THREADING C++11
   CACHE_SIZE auto
   VECTORIZATION on
)
