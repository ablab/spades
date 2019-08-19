############################################################################
# Copyright (c) 2019 University of Warwick
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(bin-analyzer CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(bin-analyzer
               main.cpp)

target_link_libraries(bin-analyzer graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(bin-analyzer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS bin-analyzer
        DESTINATION bin
        COMPONENT runtime)
