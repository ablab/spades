############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(paired_info CXX)

add_library(paired_info STATIC
            distance_estimation.cpp weighted_distance_estimation.cpp smoothing_distance_estimation.cpp
            paired_info_utils.cpp distance_estimation_utils.cpp)
target_link_libraries(paired_info modules)
