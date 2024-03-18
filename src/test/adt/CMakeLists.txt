############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2018-2022 Saint Petersburg State University
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

add_executable(phm_test
               phm_test.cpp)
target_link_libraries(phm_test utils ${COMMON_LIBRARIES} gtest)
