############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(stages-mpi CXX)

set(STAGES_MPI_SRC construction_mpi.cpp test_mpi.cpp)

add_library(stages-mpi STATIC
            ${STAGES_MPI_SRC})
