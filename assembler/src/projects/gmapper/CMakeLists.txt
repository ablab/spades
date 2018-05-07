############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(gmapper CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(spades-gmapper
               main.cpp)

target_link_libraries(spades-gmapper spades-stages graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-gmapper PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-gmapper
        DESTINATION bin
        COMPONENT runtime)
