############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(edge-profiles CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(spades-edge-profiles
               main.cpp)

target_link_libraries(spades-edge-profiles graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-edge-profiles PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

