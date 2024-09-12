############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(utils_mpi CXX)

set(utils_mpi_src logger/mpi_log_writers.cpp)

add_library(utils_mpi STATIC ${utils_mpi_src})
target_link_libraries(utils_mpi ${COMMON_LIBRARIES})
