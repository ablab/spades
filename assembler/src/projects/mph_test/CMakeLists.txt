############################################################################
# Copyright (c) 2016 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades-kmercount CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(spades-kmercount
               main.cpp)

target_link_libraries(spades-kmercount common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-kmercount PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()
