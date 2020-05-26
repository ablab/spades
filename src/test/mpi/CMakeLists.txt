############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(mpi_test CXX)

add_executable(mpi_test mpi_test.cpp)
target_link_libraries(mpi_test utils ${COMMON_LIBRARIES})
if (MPI_FOUND)
  target_link_libraries(mpi_test ${MPI_LIBRARIES})
endif()


