############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades CXX)

add_executable(spades
 main.cpp)

target_link_libraries(spades ${COMMON_LIBRARIES} input mph_index)

if (SPADES_STATIC_BUILD)
  set_target_properties(spades PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades
        RUNTIME DESTINATION bin)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../configs/debruijn"
        DESTINATION configs
        FILES_MATCHING PATTERN "*.info.template")
