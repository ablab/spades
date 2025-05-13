############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(mpi_test CXX)

set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
find_package(MPI)
if (MPI_FOUND)
  # Determine MPI vendor and MPI runtime version
  # configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/MPIVendorName.c.in"
  #                "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/MPIVendorName.c"
  #                IMMEDIATE @ONLY)
  # try_run(MPI_VENDOR_NAME_RUN MPI_HAVE_VENDOR_NAME
  #         ${CMAKE_BINARY_DIR}
  #         "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/MPIVendorName.c"
  #         RUN_OUTPUT_VARIABLE MPI_RUNTIME_NAME)
  # configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/MPIVendorVersion.c.in"
  #                "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/MPIVendorVersion.c"
  #                IMMEDIATE @ONLY)
  # try_run(MPI_VENDOR_VERSION_RUN MPI_HAVE_VENDOR_VERSION
  #         ${CMAKE_BINARY_DIR}
  #         "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/MPIVendorVersion.c"
  #         RUN_OUTPUT_VARIABLE MPI_RUNTIME_VERSION)
  message(STATUS "Detected MPI runtime: ${MPI_C_LIBRARY_VERSION_STRING}")

  if ("${MPI_C_LIBRARY_VERSION_STRING}" MATCHES "^Open MPI")
    string(REGEX REPLACE "Open MPI v([0-9]+).*" "\\1" OPENMPI_MAJOR_VERSION "${MPI_C_LIBRARY_VERSION_STRING}")
    message(STATUS "Open MPI runtime detected, major version: ${OPENMPI_MAJOR_VERSION}")
    if (OPENMPI_MAJOR_VERSION STREQUAL 3)
      message(FATAL_ERROR "Open MPI version ${OPENMPI_MAJOR_VERSION}.x is known to be buggy")
    endif()
  endif()
endif()

if (NOT MPI_ENABLE AND (NOT MPI_FOUND OR MPI_DISABLE))
  message(FATAL_ERROR "hpcSPAdees requires MPI to be enabled")
endif()

add_executable(mpi_test mpi_test.cpp)
include_directories(${SPADES_MAIN_PROJ_DIR}/hpcspades)
target_link_libraries(mpi_test utils_mpi utils llvm-support ${COMMON_LIBRARIES} ${MPI_LIBRARIES})

