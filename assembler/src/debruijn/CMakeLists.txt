############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades CXX)

add_executable(spades
               kmer_coverage_model.cpp
               config_struct.cpp
               main.cpp)

target_link_libraries(spades ${COMMON_LIBRARIES} input mph_index otkpp)

if (SPADES_STATIC_BUILD)
  set_target_properties(spades PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/debruijn"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
