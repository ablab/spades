# -*- cmake -*-

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  # Require at least gcc 9.1
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    message(FATAL_ERROR "SPAdes requires gcc version 9.1 or later")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
    message(FATAL_ERROR "SPAdes requires clang version 5 or later")
  endif()
else()
  message(WARNING "Unsupported compiler is detected. SPAdes compilation was not tested on it and may fail")
endif()

find_package(OpenMP COMPONENTS CXX)
find_package(BZip2 REQUIRED)
find_package(Readline QUIET)
set(CURSES_NEED_NCURSES TRUE)
find_package(Curses QUIET)

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

# Use included boost unless explicitly specified
if (NOT SPADES_BOOST_ROOT)
  set(BOOST_ROOT "${EXT_DIR}/include")
else()
  set(BOOST_ROOT SPADES_BOOST_ROOT)
endif()
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED)
