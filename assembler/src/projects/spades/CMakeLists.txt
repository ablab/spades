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
            chromosome_removal.cpp
            extract_domains.cpp
            domain_matcher.cpp
            domain_graph_construction.cpp
            restricted_edges_filling.cpp
            load_graph.cpp
            domain_graph.cpp)

target_link_libraries(spades-stages hmmercpp spoa easel)

add_executable(spades-core
               pipeline.cpp
               main.cpp
               series_analysis.cpp
               ../mts/contig_abundance.cpp)

target_link_libraries(spades-core spades-stages graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-core PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-core
        DESTINATION bin
        COMPONENT spades)
# Configs
install(DIRECTORY "configs/" # Trailing / is important
        DESTINATION share/spades/configs/debruijn
        COMPONENT spades
        FILES_MATCHING PATTERN "*.info")
# biosyntheticSPAdes default HMM profiles
install(DIRECTORY "biosynthetic_spades_hmms"
        DESTINATION share/spades
        COMPONENT spades)
# coronaSPAdes default HMM profiles
install(DIRECTORY "coronaspades_hmms"
        DESTINATION share/spades
        COMPONENT spades)

# Main pipeline script
add_subdirectory(pipeline)

# Test datasets
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test_dataset"
        DESTINATION share/spades
        COMPONENT spades-test)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test_dataset_plasmid"
        DESTINATION share/spades
        COMPONENT spades-test)
