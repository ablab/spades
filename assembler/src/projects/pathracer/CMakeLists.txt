############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(pathracer CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(pathracer-core STATIC
            debruijn_graph_cursor.cpp fees.cpp
            find_best_path.cpp
            fasta_reader.cpp)
target_link_libraries(pathracer-core hmmercpp assembly_graph common_modules)

add_executable(pathracer
               main.cpp pathracer.cpp)
target_link_libraries(pathracer
                      pathracer-core
                      graphio utils ${COMMON_LIBRARIES})
# target_compile_options(pathracer PRIVATE -g)
install(TARGETS pathracer
        DESTINATION bin
        COMPONENT runtime)

if (SPADES_STATIC_BUILD)
  set_target_properties(pathracer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

add_executable(align_fs
               main.cpp align_fs.cpp)
target_link_libraries(align_fs
                      pathracer-core
                      graphio utils ${COMMON_LIBRARIES})

add_executable(pathracer_align_kmers_test
               main.cpp align_kmers.cpp)
target_link_libraries(pathracer_align_kmers_test
                      pathracer-core
                      graphio utils ${COMMON_LIBRARIES})
set_target_properties(pathracer_align_kmers_test PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

add_executable(graph_exact_finder graph_exact_finder.cpp)
target_link_libraries(graph_exact_finder
                      pathracer-core
                      graphio utils ${COMMON_LIBRARIES})
set_target_properties(graph_exact_finder PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

add_executable(open_event_graph open_event_graph.cpp)
target_link_libraries(open_event_graph
                      pathracer-core
                      graphio utils ${COMMON_LIBRARIES})
set_target_properties(open_event_graph PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)


add_library(gtest_main_segfault_handler gtest_main.cpp)
target_link_libraries(gtest_main_segfault_handler gtest input utils ${COMMON_LIBRARIES})

add_executable(pathracer-test-levenshtein find_best_path.cpp fees.cpp test.cpp)
target_link_libraries(pathracer-test-levenshtein gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})
add_test(NAME pathracer-levenshtein COMMAND pathracer-test-levenshtein)

add_executable(pathracer-test-aa test-aa.cpp graph.cpp fees.cpp)
target_link_libraries(pathracer-test-aa gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})
add_test(NAME pathracer-aa COMMAND pathracer-test-aa)

add_executable(pathracer-test-depth-int test-depth.cpp graph.cpp fees.cpp)
target_link_libraries(pathracer-test-depth-int gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})
add_test(NAME pathracer-depth-int COMMAND pathracer-test-depth-int)
add_executable(pathracer-test-cursor-utils test-cursor-utils.cpp graph.cpp fees.cpp)
target_link_libraries(pathracer-test-cursor-utils gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})
add_test(NAME pathracer-cursor-utils COMMAND pathracer-test-cursor-utils)
# add_executable(pathracer-test-stack-limit test-stack-limit.cpp graph.cpp fees.cpp)
# target_link_libraries(pathracer-test-stack-limit gtest_main_segfault_handler hmmercpp input utils pipeline ${COMMON_LIBRARIES})
# add_test(NAME pathracer-stack-limit COMMAND pathracer-test-stack-limit)
