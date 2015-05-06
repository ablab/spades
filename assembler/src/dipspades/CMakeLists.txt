############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(dipspades CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/debruijn)

add_executable(dipspades
               ../debruijn/kmer_coverage_model.cpp
               ../debruijn/config_struct.cpp
               ../debruijn/path_extend/pe_config_struct.cpp
               ../debruijn/stage.cpp
               ../debruijn/construction.cpp
               ../debruijn/gap_closer.cpp
               ../debruijn/simplification.cpp
               ../debruijn/mismatch_correction.cpp
               ../debruijn/pair_info_count.cpp
               ../debruijn/distance_estimation.cpp
               ../debruijn/repeat_resolving.cpp
               ../debruijn/genomic_info_filler.cpp
		  dipspades_config.cpp
	      utils/files_utils.cpp
          main.cpp)

target_link_libraries(dipspades input mph_index BamTools yaml-cpp nlopt ${COMMON_LIBRARIES} ssw)

if (SPADES_STATIC_BUILD)
  set_target_properties(dipspades PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS dipspades
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/dipspades"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
