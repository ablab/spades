############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# Copyright (c) 2019 University of Warwick
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(unitig-coverage CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(unitig-coverage
               main.cpp profile_storage.cpp)

target_link_libraries(unitig-coverage graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(unitig-coverage PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

