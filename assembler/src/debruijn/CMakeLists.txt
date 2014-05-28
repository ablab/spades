############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades CXX)

add_executable(spades
               kmer_coverage_model.cpp
               config_struct.cpp
               path_extend/pe_config_struct.cpp
               stage.cpp
               construction.cpp
               gap_closer.cpp
               simplification.cpp
               mismatch_correction.cpp
               pair_info_count.cpp
               distance_estimation.cpp
               repeat_resolving.cpp
               genomic_info_filler.cpp
	           pacbio_aligning.cpp
               main.cpp)

# We have to do this to minimize changes in ConsensusCore itself
include_directories(${EXT_DIR}/include/ConsensusCore)

target_link_libraries(spades ConsensusCore input mph_index nlopt BamTools yaml-cpp ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/debruijn"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
