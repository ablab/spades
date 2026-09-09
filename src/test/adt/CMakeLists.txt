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

add_executable(adt_test
               loser_tree_test.cpp phm_test.cpp
               test.cpp)
target_link_libraries(adt_test utils ${COMMON_LIBRARIES} teamcity_gtest gtest)

add_test(NAME adt_test COMMAND adt_test)
