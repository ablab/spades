############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades CXX)

add_library(debruijn STATIC
            ../include/sam/read.cpp
            ../include/sam/sam_reader.cpp
            kmer_coverage_model.cpp
            config_struct.cpp
            path_extend/pe_config_struct.cpp
            path_extend/bidirectional_path.cpp
            path_extend/scaff_supplementary.cpp
            path_extend/extension_chooser2015.cpp
            genome_consistance_checker.cpp
            stage.cpp
            construction.cpp
            gap_closer.cpp
            simplification.cpp
            mismatch_correction.cpp
            pair_info_count.cpp
            second_phase_setup.cpp
            distance_estimation.cpp
            repeat_resolving.cpp
            genomic_info_filler.cpp
            pacbio_aligning.cpp
            genome_storage.cpp 
            path_extend/scaffold_graph.cpp 
            path_extend/scaffold_graph_visualizer.cpp 
            path_extend/connection_condition2015.cpp
            bwa_pair_info_count.cpp
            bwa_pair_info_filler.cpp)
  
# We have to do this to minimize changes in ConsensusCore itself
target_include_directories(debruijn PRIVATE ${EXT_DIR}/include/ConsensusCore)
target_link_libraries(debruijn ConsensusCore input samtools mph_index nlopt BamTools yaml-cpp ${COMMON_LIBRARIES})

add_executable(spades
               main.cpp)
target_link_libraries(spades debruijn)

if (SPADES_STATIC_BUILD)
  set_target_properties(spades PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/debruijn"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
