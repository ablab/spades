############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(contig_corrector CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(contig_corrector
               contig_replacer/contig_replacer.cpp
               error_analyzer/error_analyzer.cpp
               error_analyzer/sequence_fragments.cpp
               error_analyzer/statistics.cpp
               helpers/parser.cpp
               helpers/replacer.cpp
               main.cpp
               sequence_clusterer/sequence_clusterer.cpp
               sequence_corrector/filler_chooser.cpp
               sequence_corrector/launcher.cpp
               sequence_corrector/sequence_corrector.cpp
              )

target_link_libraries(contig_corrector graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(contig_corrector PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS contig_corrector
        DESTINATION bin
        COMPONENT runtime)
