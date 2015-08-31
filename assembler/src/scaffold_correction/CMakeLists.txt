############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(moleculo CXX)


add_executable(scaffold_correction
               main.cpp)
target_include_directories(scaffold_correction PRIVATE ${CMAKE_SOURCE_DIR}/debruijn)
target_link_libraries(scaffold_correction debruijn)

if (SPADES_STATIC_BUILD)
  set_target_properties(scaffold_correction PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS scaffold_correction
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/scaffold_correction"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
