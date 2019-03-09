############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(pathracer CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(pathracer
               omnigraph_wrapper.cpp fees.cpp
               find_best_path.cpp
               main.cpp
               fasta_reader.cpp
               pathracer.cpp)
target_link_libraries(pathracer graphio modules hmmercpp input utils pipeline ${COMMON_LIBRARIES} hattrie cityhash modules assembly_graph common_modules)
# target_compile_options(pathracer PRIVATE -g)

add_executable(pathracer_align_kmers_test
               fasta_reader.cpp
               omnigraph_wrapper.cpp fees.cpp
               find_best_path.cpp
               main.cpp
               align_kmers.cpp)
target_link_libraries(pathracer_align_kmers_test graphio modules hmmercpp input utils pipeline ${COMMON_LIBRARIES} hattrie cityhash modules assembly_graph common_modules)
# set_target_properties(pathracer_align_kmers_test PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

add_executable(graph_exact_finder
               fasta_reader.cpp
               omnigraph_wrapper.cpp fees.cpp
               find_best_path.cpp
               graph_exact_finder.cpp)
target_link_libraries(graph_exact_finder graphio modules hmmercpp input utils pipeline ${COMMON_LIBRARIES} hattrie cityhash modules assembly_graph common_modules)

add_executable(open_event_graph
               fasta_reader.cpp
               omnigraph_wrapper.cpp fees.cpp
               find_best_path.cpp
               open_event_graph.cpp)
target_link_libraries(open_event_graph graphio modules hmmercpp input utils pipeline ${COMMON_LIBRARIES} hattrie cityhash modules assembly_graph common_modules)

install(TARGETS pathracer
        DESTINATION bin
        COMPONENT runtime)

if (SPADES_STATIC_BUILD)
  set_target_properties(pathracer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

add_library(gtest_main_segfault_handler gtest_main.cpp)
target_link_libraries(gtest_main_segfault_handler gtest input utils ${COMMON_LIBRARIES})

add_executable(graph-hmm-test-levenshtein graph.cpp fees.cpp test.cpp)
target_link_libraries(graph-hmm-test-levenshtein gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})
add_test(NAME pathracer-levenshtein COMMAND graph-hmm-test-levenshtein)

add_executable(graph-hmm-test-aa test-aa.cpp graph.cpp fees.cpp)
target_link_libraries(graph-hmm-test-aa gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})
add_test(NAME pathracer-aa COMMAND graph-hmm-test-aa)

# add_executable(graph-hmm-test-stack-limit test-stack-limit.cpp graph.cpp fees.cpp)
# target_link_libraries(graph-hmm-test-stack-limit gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})
# add_test(NAME pathracer-stack-limit COMMAND graph-hmm-test-stack-limit)
