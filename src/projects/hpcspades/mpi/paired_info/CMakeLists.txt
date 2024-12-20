############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(paired_info_mpi CXX)

add_library(paired_info_mpi STATIC
            distance_estimation.cpp
            distance_estimation_utils.cpp)

target_link_libraries(paired_info_mpi modules)
