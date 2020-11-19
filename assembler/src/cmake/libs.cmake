# -*- cmake -*-

# Collect all the necessary common libraries
set(COMMON_LIBRARIES ${SYSTEM_LIBRARIES})

if (SPADES_USE_JEMALLOC)
  set(COMMON_LIBRARIES jemalloc-static ${COMMON_LIBRARIES})
endif()

# Add TCMalloc
if (SPADES_USE_TCMALLOC)
  set(COMMON_LIBRARIES ${TCMALLOC_LIBRARIES} ${COMMON_LIBRARIES})
endif()

# Add format library
set(COMMON_LIBRARIES format ${COMMON_LIBRARIES})
# Add version
set(COMMON_LIBRARIES version ${COMMON_LIBRARIES})

# Really, should be last here
if (SPADES_USE_MIMALLOC)
  set(COMMON_LIBRARIES mimalloc-obj ${COMMON_LIBRARIES})
endif()
