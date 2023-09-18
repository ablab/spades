############################################################################
# Copyright (c) 2023 SPAdes authors
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# -*- cmake -*-

function(canonicalize_tool_name name output)
  string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" nameStrip ${name})
  string(REPLACE "-" "_" nameUNDERSCORE ${nameStrip})
  string(TOUPPER ${nameUNDERSCORE} nameUPPER)
  set(${output} "${nameUPPER}" PARENT_SCOPE)
endfunction(canonicalize_tool_name)

# Custom add_subdirectory wrapper
# Takes in a project name (i.e. SPAdes), the subdirectory name, and an optional
# path if it differs from the name.
function(add_spades_subdirectory project type name)
  set(add_spades_external_dir "${ARGN}")
  if ("${add_spades_external_dir}" STREQUAL "")
    set(add_spades_external_dir ${name})
  endif()
  canonicalize_tool_name(${name} nameUPPER)
  set(canonical_full_name ${project}_${type}_${nameUPPER})
  get_property(already_processed GLOBAL PROPERTY ${canonical_full_name}_PROCESSED)
  if (already_processed)
    return()
  endif()
  set_property(GLOBAL PROPERTY ${canonical_full_name}_PROCESSED YES)

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${add_spades_external_dir}/CMakeLists.txt)
    # Treat it as in-tree subproject.
    option(${canonical_full_name}_BUILD
           "Whether to build ${name} as part of ${project}" On)
    mark_as_advanced(${project}_${type}_${name}_BUILD)
    if (${canonical_full_name}_BUILD)
      add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${add_spades_external_dir} ${add_spades_external_dir})
    endif()
  else()
    set(SPADES_EXTERNAL_${nameUPPER}_SOURCE_DIR
      "${SPADES_EXTERNAL_${nameUPPER}_SOURCE_DIR}"
      CACHE PATH "Path to ${name} source directory")
    set(${canonical_full_name}_BUILD_DEFAULT ON)
    if (NOT SPADES_EXTERNAL_${nameUPPER}_SOURCE_DIR OR NOT EXISTS ${SPADES_EXTERNAL_${nameUPPER}_SOURCE_DIR})
      set(${canonical_full_name}_BUILD_DEFAULT OFF)
    endif()
    if ("${SPADES_EXTERNAL_${nameUPPER}_BUILD}" STREQUAL "OFF")
      set(${canonical_full_name}_BUILD_DEFAULT OFF)
    endif()
    option(${canonical_full_name}_BUILD
      "Whether to build ${name} as part of SPAdes"
      ${${canonical_full_name}_BUILD_DEFAULT})
    if (${canonical_full_name}_BUILD)
      if (EXISTS ${SPADES_EXTERNAL_${nameUPPER}_SOURCE_DIR})
        add_subdirectory(${SPADES_EXTERNAL_${nameUPPER}_SOURCE_DIR} ${add_spades_external_dir})
      elseif(NOT "${SPADES_EXTERNAL_${nameUPPER}_SOURCE_DIR}" STREQUAL "")
        message(WARNING "Nonexistent directory for ${name}: ${SPADES_EXTERNAL_${nameUPPER}_SOURCE_DIR}")
      endif()
    endif()
  endif()
endfunction()

# Add external project that may want to be built as part of SPAdes such as SPAligner,
# Pathracer, and BinSPreader. This adds two options. One for the source directory of the
# project, which defaults to ${CMAKE_CURRENT_SOURCE_DIR}/${name}. Another to
# enable or disable building it with everything else.
# Additional parameter can be specified as the name of directory.
macro(add_spades_external_project name)
  add_spades_subdirectory(SPADES TOOL ${name} ${ARGN})
endmacro()

macro(add_spades_tool_subdirectory name)
  add_spades_external_project(${name})
endmacro(add_spades_tool_subdirectory)

function(create_subdirectory_options project type)
  file(GLOB sub-dirs "${CMAKE_CURRENT_SOURCE_DIR}/*")
  foreach(dir ${sub-dirs})
    if (IS_DIRECTORY "${dir}" AND EXISTS "${dir}/CMakeLists.txt")
      canonicalize_tool_name(${dir} name)
      option(${project}_${type}_${name}_BUILD
           "Whether to build ${name} as part of ${project}" On)
      mark_as_advanced(${project}_${type}_${name}_BUILD)
    endif()
  endforeach()
endfunction(create_subdirectory_options)

function(create_spades_tool_options)
  create_subdirectory_options(SPADES TOOL)
endfunction(create_spades_tool_options)

function(spades_add_implicit_projects project)
  set(list_of_implicit_subdirs "")
  file(GLOB sub-dirs "${CMAKE_CURRENT_SOURCE_DIR}/*")
  foreach(dir ${sub-dirs})
    if(IS_DIRECTORY "${dir}" AND EXISTS "${dir}/CMakeLists.txt")
      canonicalize_tool_name(${dir} name)
      # I don't like special casing things by order, but the spades-driver ends up
      # linking the object libraries from all the tools that opt-in, so adding
      # it separately at the end is probably the simplest case.
      if("${name}" STREQUAL "SPADES_DRIVER")
        continue()
      endif()
      if (${project}_TOOL_${name}_BUILD)
        get_filename_component(fn "${dir}" NAME)
        list(APPEND list_of_implicit_subdirs "${fn}")
      endif()
    endif()
  endforeach()

  foreach(external_proj ${list_of_implicit_subdirs})
    add_spades_subdirectory(${project} TOOL "${external_proj}" ${ARGN})
  endforeach()
endfunction(spades_add_implicit_projects)

function(add_spades_implicit_projects)
  spades_add_implicit_projects(SPADES)
endfunction(add_spades_implicit_projects)
