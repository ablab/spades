############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2019-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(debruijn_test CXX)

add_executable(debruijn_test
               graph_core_test.cpp histogram_test.cpp paired_info_test.cpp overlap_analysis_test.cpp
               simplification_test.cpp test_utils.cpp construction_test.cpp io_test.cpp
               path_extend_test.cpp graphio.cpp overlap_removal_test.cpp graph_alignment_test.cpp v_overlaps.cpp
               test.cpp)
target_link_libraries(debruijn_test common_modules input ${COMMON_LIBRARIES} graphio teamcity_gtest gtest)
add_test(NAME debruijn_test COMMAND debruijn_test)
