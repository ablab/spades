############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(graph-hmm CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(graph-hmm
               hmmfile.cpp hmmmatcher.cpp
               main.cpp)

add_executable(looper
               looper.cpp
               hmmfile.cpp
               demo.cpp
               graph.cpp)

target_link_libraries(graph-hmm easel hmmer input utils pipeline ${COMMON_LIBRARIES})
target_link_libraries(looper easel hmmer input utils pipeline ${COMMON_LIBRARIES})

add_subdirectory(googletest)
add_library(gtest_main_segfault_handler gtest_main.cpp)
target_link_libraries(gtest_main_segfault_handler gtest input utils ${COMMON_LIBRARIES})

add_executable(looper-test1 graph.cpp demo.cpp hmmfile.cpp demo.cpp test.cpp)
target_link_libraries(looper-test1 gtest_main_segfault_handler easel hmmer input utils pipeline ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(graph-hmm PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(looper PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()
