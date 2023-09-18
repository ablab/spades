############################################################################
# Copyright (c) 2023 SPAdes authors
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades-tools CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(spades-convert-bin-to-fasta
               convert_bin_to_fasta.cpp)
add_executable(spades-gbuilder
               gbuilder.cpp)
add_executable(spades-gmapper
               gmapper.cpp)
add_executable(spades-gsimplifier
               ../cds_subgraphs/profile_storage.cpp
               gsimplifier.cpp)
add_executable(spades-kmercount
               kmercount.cpp)
add_executable(spades-read-filter
               read_filter.cpp)
add_executable(spades-kmer-estimating
               kmer_estimating.cpp)

target_link_libraries(spades-convert-bin-to-fasta common_modules ${COMMON_LIBRARIES})
target_link_libraries(spades-gbuilder graphio common_modules ${COMMON_LIBRARIES})
target_link_libraries(spades-gmapper spades-stages graphio common_modules ${COMMON_LIBRARIES})
target_link_libraries(spades-gsimplifier toolchain common_modules ${COMMON_LIBRARIES})
target_link_libraries(spades-kmercount common_modules ${COMMON_LIBRARIES})
target_link_libraries(spades-read-filter common_modules ${COMMON_LIBRARIES})
target_link_libraries(spades-kmer-estimating common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-convert-bin-to-fasta PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(spades-gbuilder PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(spades-gmapper PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(spades-gsimplifier PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(spades-kmercount PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(spades-read-filter PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-convert-bin-to-fasta
        DESTINATION bin
        COMPONENT tools)

install(TARGETS spades-gbuilder
        DESTINATION bin
        COMPONENT tools)

install(TARGETS spades-gmapper
        DESTINATION bin
        COMPONENT tools)      

install(TARGETS spades-gsimplifier
        DESTINATION bin
        COMPONENT tools)      

install(TARGETS spades-kmercount
        DESTINATION bin
        COMPONENT tools)

install(TARGETS spades-read-filter
        DESTINATION bin
        COMPONENT tools)

install(TARGETS spades-kmer-estimating
        DESTINATION bin
        COMPONENT tools)
