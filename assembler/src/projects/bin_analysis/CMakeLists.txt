############################################################################
# Copyright (c) 2019 University of Warwick
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(bin-analysis CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(mag-improve
               main.cpp)

target_link_libraries(mag-improve graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(mag-improve PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS mag-improve
        DESTINATION bin
        COMPONENT runtime)
