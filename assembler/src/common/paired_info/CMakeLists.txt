############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(paired_info CXX)

add_library(paired_info STATIC
        distance_estimation.cpp weighted_distance_estimation.cpp smoothing_distance_estimation.cpp)
