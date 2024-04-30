############################################################################
# Copyright (c) 2023 SPAdes authors
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# -*- cmake -*-

get_cmake_property(vars CACHE_VARIABLES)
foreach(var ${vars})
  get_property(var_type CACHE ${var} PROPERTY TYPE)
  if (NOT var_type STREQUAL "STATIC" AND
      NOT var_type STREQUAL "INTERNAL" AND
      NOT var STREQUAL "CMAKE_INSTALL_PREFIX")
      list(APPEND CL_ARGS "-D${var}:${var_type}=${${var}}")
    endif()
endforeach()

include(ExternalProject)

ExternalProject_Add(ncbi_vdb_ext
                    GIT_REPOSITORY "https://github.com/ncbi/ncbi-vdb.git"
                    GIT_TAG "3.1.0"
                    GIT_SHALLOW 1
                    PATCH_COMMAND patch -p1 < ${CMAKE_CURRENT_LIST_DIR}/vdb.patch
                    BUILD_ALWAYS OFF
                    UPDATE_DISCONNECTED 1
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                    CMAKE_CACHE_ARGS ${CL_ARGS})
ExternalProject_Get_property(ncbi_vdb_ext INSTALL_DIR)
set(NCBI_VDB_PREFIX  ${INSTALL_DIR})
set(NCBI_VDB_INCLUDE ${NCBI_VDB_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
set(NCBI_VDB_LIB     ${NCBI_VDB_PREFIX}/${CMAKE_INSTALL_LIBDIR})

ExternalProject_Add(ncbi_sratools_ext
                    GIT_REPOSITORY "https://github.com/ncbi/sra-tools.git"
                    GIT_TAG "3.1.0"
                    GIT_SHALLOW 1
                    PATCH_COMMAND patch -p1 < ${CMAKE_CURRENT_LIST_DIR}/sratools.patch
                    BUILD_ALWAYS OFF
                    UPDATE_DISCONNECTED 1
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                               -DVDB_INCDIR=${NCBI_VDB_INCLUDE}
                               -DVDB_LIBDIR=${NCBI_VDB_LIB}
                               -DNO_JAVA=ON
                    CMAKE_CACHE_ARGS ${CL_ARGS})
add_dependencies(ncbi_sratools_ext ncbi_vdb_ext)
ExternalProject_Get_property(ncbi_sratools_ext INSTALL_DIR)
set(NCBI_SRATOOLS_PREFIX  ${INSTALL_DIR})
set(NCBI_SRATOOLS_INCLUDE ${NCBI_SRATOOLS_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
set(NCBI_SRATOOLS_LIB     ${NCBI_SRATOOLS_PREFIX}/${CMAKE_INSTALL_LIBDIR})

add_library(ncbi-sdk INTERFACE)

add_library(ncbi-vdb STATIC IMPORTED)
set_target_properties(ncbi-vdb PROPERTIES
                      IMPORTED_LOCATION ${NCBI_VDB_LIB}/${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-vdb-static${CMAKE_STATIC_LIBRARY_SUFFIX}
                      INTERFACE_INCLUDE_DIRECTORIES ${NCBI_VDB_INCLUDE})
add_dependencies(ncbi-vdb ncbi_vdb_ext)

add_library(ncbi-c++ STATIC IMPORTED GLOBAL)
set_target_properties(ncbi-c++ PROPERTIES
                      IMPORTED_LOCATION ${NCBI_SRATOOLS_LIB}/${CMAKE_STATIC_LIBRARY_PREFIX}ngs-c++-static${CMAKE_STATIC_LIBRARY_SUFFIX}
                      INTERFACE_INCLUDE_DIRECTORIES ${NCBI_SRATOOLS_INCLUDE})
add_dependencies(ncbi-c++ ncbi_sratools_ext)

add_library(ncbi-ngs STATIC IMPORTED)
set_target_properties(ncbi-ngs PROPERTIES
                      IMPORTED_LOCATION ${NCBI_SRATOOLS_LIB}/${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-ngs-static${CMAKE_STATIC_LIBRARY_SUFFIX}
                      INTERFACE_INCLUDE_DIRECTORIES ${NCBI_SRATOOLS_INCLUDE})
add_dependencies(ncbi-ngs ncbi_sratools_ext)

add_library(ncbi-ngs-c++ STATIC IMPORTED)
set_target_properties(ncbi-ngs-c++ PROPERTIES
                      IMPORTED_LOCATION ${NCBI_SRATOOLS_LIB}/${CMAKE_STATIC_LIBRARY_PREFIX}ncbi-ngs-c++-static${CMAKE_STATIC_LIBRARY_SUFFIX}
                      INTERFACE_INCLUDE_DIRECTORIES ${NCBI_SRATOOLS_INCLUDE})
add_dependencies(ncbi-ngs-c++ ncbi_statools_ext)


set(NCBI_SDK_INCLUDES ${NCBI_SRATOOLS_INCLUDE} ${NCBI_VDB_INCLUDE})
set(NCBI_SDK_LIBS ncbi-ngs-c++ ncbi-ngs ncbi-vdb ncbi-c++)

set_target_properties(ncbi-sdk PROPERTIES
                      INTERFACE_LINK_LIBRARIES "${NCBI_SDK_LIBS}"
                      INTERFACE_INCLUDE_DIRECTORIES "${NCBI_SDK_INCLUDES}")
add_dependencies(ncbi-sdk ncbi-c++ ncbi-vdb ncbi-ngs ncbi-ngs-с++)
