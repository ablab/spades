############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades CXX)

add_library(deb STATIC
            gap_closer.cpp
            simplification.cpp
            mismatch_correction.cpp
            pair_info_count.cpp
            second_phase_setup.cpp
            distance_estimation.cpp
            repeat_resolving.cpp
            pacbio_aligning.cpp
            chromosome_removal.cpp)
  
target_include_directories(deb PRIVATE ${EXT_DIR}/include/ConsensusCore)
target_link_libraries(deb ConsensusCore spades_modules cityhash nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})

add_executable(spades main.cpp)
target_link_libraries(spades deb)

if (SPADES_STATIC_BUILD)
  set_target_properties(spades PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/debruijn"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
