############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(ionhammer CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ionhammer
               kmer_data.cpp
               hamcluster.cpp
               subcluster.cpp
               err_helper_table.cpp
               config_struct.cpp
               expander.cpp
               seqeval/BaseHypothesisEvaluator.cpp
               seqeval/TreephaserLite.cpp
               main.cpp)

target_link_libraries(ionhammer input utils pipeline mph_index BamTools ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(ionhammer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS ionhammer
        RUNTIME DESTINATION bin)
install(DIRECTORY "${SPADES_CFG_DIR}/ionhammer"
        DESTINATION share/spades/configs
                FILES_MATCHING PATTERN "*.cfg")
