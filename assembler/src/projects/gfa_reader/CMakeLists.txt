############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(gfa-reader CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(spades-gfareader
               main.cpp)

target_link_libraries(spades-gfareader gfaio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-gfareader PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

