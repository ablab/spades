############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades-bin-converter CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(spades-convert-bin-to-fasta
               convert_bin_to_fasta.cpp)

target_link_libraries(spades-convert-bin-to-fasta common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-convert-bin-to-fasta PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-convert-bin-to-fasta
        DESTINATION bin
        COMPONENT runtime)