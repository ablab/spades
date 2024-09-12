############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(alignment_mpi CXX)

add_library(alignment_mpi STATIC
            sequence_mapper_notifier_mpi.cpp)

target_link_libraries(alignment_mpi modules)
