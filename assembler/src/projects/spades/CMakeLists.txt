############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades CXX)

add_executable(spades main.cpp
            gap_closer.cpp
            mismatch_correction.cpp
            pair_info_count.cpp
            second_phase_setup.cpp
            distance_estimation.cpp
            repeat_resolving.cpp
            contig_output_stage.cpp
            hybrid_aligning.cpp
            chromosome_removal.cpp
            ../mts/contig_abundance.cpp)

target_include_directories(spades PRIVATE ${EXT_DIR}/include/ConsensusCore)
target_link_libraries(spades ConsensusCore common_modules nlopt BamTools ssw ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/debruijn"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
