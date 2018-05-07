############################################################################
# Copyright (c) 2017 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades-gbuilder CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(spades-gbuilder
               main.cpp)

target_link_libraries(spades-gbuilder graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-gbuilder PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-gbuilder
        DESTINATION bin
        COMPONENT runtime)
