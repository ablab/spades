############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(pathracer CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(pathracer
               omnigraph_wrapper.cpp fees.cpp
               main.cpp)

add_executable(looper
               looper.cpp
               fasta_reader.cpp fees.cpp
               graph.cpp)
set_target_properties(looper PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

target_link_libraries(pathracer graphio modules hmmercpp input utils pipeline ${COMMON_LIBRARIES})
target_link_libraries(looper hmmercpp input utils pipeline ${COMMON_LIBRARIES})

add_subdirectory(googletest)
add_library(gtest_main_segfault_handler gtest_main.cpp)
target_link_libraries(gtest_main_segfault_handler gtest input utils ${COMMON_LIBRARIES})

add_executable(graph-hmm-test-levenshtein graph.cpp fees.cpp test.cpp)
target_link_libraries(graph-hmm-test-levenshtein gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})

add_executable(graph-hmm-test-aa test-aa.cpp graph.cpp fees.cpp)
target_link_libraries(graph-hmm-test-aa gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})

add_executable(graph-hmm-test-stack-limit test-stack-limit.cpp graph.cpp fees.cpp)
target_link_libraries(graph-hmm-test-stack-limit gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(pathracer PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(looper PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()
