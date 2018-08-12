############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(adt_test CXX)

add_executable(pack_demo
               pack_demo.cpp)
target_link_libraries(pack_demo ${COMMON_LIBRARIES})

add_executable(loser_tree_test
               loser_tree_test.cpp)
target_link_libraries(loser_tree_test ${COMMON_LIBRARIES})
