############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(pipeline_mpi CXX)

add_library(pipeline_mpi STATIC mpi_stage.cpp)

target_link_libraries(pipeline_mpi pipeline ${MPI_LIBRARIES})
