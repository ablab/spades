# -*- cmake -*-

# Collect all the necessary common libraries
set(COMMON_LIBRARIES ${ZLIB_LIBRARIES} ${SYSTEM_LIBRARIES})

# Add TCMalloc
if (SPADES_USE_TCMALLOC)
  set(COMMON_LIBRARIES ${COMMON_LIBRARIES} ${TCMALLOC_LIBRARIES})
endif()
