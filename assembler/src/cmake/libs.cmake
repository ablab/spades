# -*- cmake -*-

# Collect all the necessary common libraries
set(COMMON_LIBRARIES  ${ZLIB_LIBRARIES} ${SYSTEM_LIBRARIES})

if (SPADES_USE_JEMALLOC)
  set(COMMON_LIBRARIES jemalloc ${COMMON_LIBRARIES})
endif()

# Add TCMalloc
if (SPADES_USE_TCMALLOC)
  set(COMMON_LIBRARIES ${TCMALLOC_LIBRARIES} ${COMMON_LIBRARIES})
endif()
