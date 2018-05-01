############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades CXX)

add_library(spades-stages STATIC
            gap_closer.cpp
            mismatch_correction.cpp
            pair_info_count.cpp
            second_phase_setup.cpp
            distance_estimation.cpp
            repeat_resolving.cpp
            contig_output_stage.cpp
            hybrid_aligning.cpp
            chromosome_removal.cpp)
target_link_libraries(spades-stages ConsensusCore)

add_executable(spades-core main.cpp
                           series_analysis.cpp
                           ../mts/contig_abundance.cpp)

target_link_libraries(spades-core spades-stages graphio common_modules nlopt BamTools ssw ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-core PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-core
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/debruijn"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
