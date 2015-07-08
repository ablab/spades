# -*- cmake -*-

# Collect all the necessary common libraries
set(COMMON_LIBRARIES ${ZLIB_LIBRARIES} ${SYSTEM_LIBRARIES})

if (SPADES_USE_JEMALLOC)
  if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(COMMON_LIBRARIES "-Wl,-force_load" jemalloc ${COMMON_LIBRARIES})
  else()
    set(COMMON_LIBRARIES jemalloc ${COMMON_LIBRARIES})
  endif()
endif()

# Add TCMalloc
if (SPADES_USE_TCMALLOC)
  set(COMMON_LIBRARIES ${TCMALLOC_LIBRARIES} ${COMMON_LIBRARIES})
endif()

# Add format library
set(COMMON_LIBRARIES format ${COMMON_LIBRARIES})
