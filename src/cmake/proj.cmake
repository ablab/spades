############################################################################
# Copyright (c) 2023 SPAdes authors
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# -*- cmake -*-

# Side-by-side subprojects layout: automatically set the
# SPADES_EXTERNAL_${project}_SOURCE_DIR using SPADES_ALL_PROJECTS
set(SPADES_ALL_PROJECTS "spades;hammer;ionhammer;corrector;spaligner;spades_tools;binspreader;pathracer;hpcspades")
set(SPADES_EXTRA_PROJECTS "mts;online_vis;cds_subgraphs")
set(SPADES_KNOWN_PROJECTS "${SPADES_ALL_PROJECTS};${SPADES_EXTRA_PROJECTS}")
set(SPADES_ENABLE_PROJECTS "" CACHE STRING
    "Semicolon-separated list of projects to build (${SPADES_KNOWN_PROJECTS}), or \"all\".")

# Make sure expansion happens first to not handle "all" in rest of the checks.
if (SPADES_ENABLE_PROJECTS STREQUAL "all")
  set(SPADES_ENABLE_PROJECTS "${SPADES_ALL_PROJECTS}")
elseif(SPADES_ENABLE_PROJECTS STREQUAL "release")
  # Exclude hpcSPAdes from release snapshots by default
  set(SPADES_RELEASE_PROJECTS "${SPADES_ALL_PROJECTS}")
  list(REMOVE_ITEM SPADES_RELEASE_PROJECTS "hpcspades")
  set(SPADES_ENABLE_PROJECTS "${SPADES_RELEASE_PROJECTS}")
endif()

# Always include SPAdes by default
if (SPADES_ENABLE_PROJECTS STREQUAL "")
  set(SPADES_ENABLE_PROJECTS "spades;spades_tools")
endif()

foreach(proj ${SPADES_ENABLE_PROJECTS})
  if (NOT "${proj}" IN_LIST SPADES_KNOWN_PROJECTS)
     MESSAGE(FATAL_ERROR "${proj} isn't a known project: ${SPADES_KNOWN_PROJECTS}.")
  endif()
endforeach()

if (SPADES_BUILD_INTERNAL)
  # Always build SPAdes for internal projects
  list(APPEND SPADES_ENABLE_PROJECTS "spades")
  list(APPEND SPADES_ENABLE_PROJECTS "mts")
  list(APPEND SPADES_ENABLE_PROJECTS "online_vis")
endif()

# Some inter-project dependencies
if ("spades" IN_LIST SPADES_ENABLE_PROJECTS)
  if (NOT "hammer" IN_LIST SPADES_ENABLE_PROJECTS)
    message(STATUS "Enabling BayesHammer as a dependency to SPAdes")
    list(APPEND SPADES_ENABLE_PROJECTS "hammer")
  endif()

  if (NOT "ionhammer" IN_LIST SPADES_ENABLE_PROJECTS)
    message(STATUS "Enabling IonHammer as a dependency to SPAdes")
    list(APPEND SPADES_ENABLE_PROJECTS "ionhammer")
  endif()

  if (NOT "corrector" IN_LIST SPADES_ENABLE_PROJECTS)
    message(STATUS "Enabling Corrector as a dependency to SPAdes")
    list(APPEND SPADES_ENABLE_PROJECTS "corrector")
  endif()
endif()

# SPADES_ENABLE_PROJECTS_USED is `ON` if the user has ever used the
# `SPADES_ENABLE_PROJECTS` CMake cache variable.  This exists for
# several reasons:
#  * As an indicator that the `SPADES_ENABLE_PROJECTS` list is now the single
#    source of truth for which projects to build. This means we will ignore user
#    supplied `SPADES_TOOL_<project>_BUILD` CMake cache variables and overwrite
#    them.
#
#  * The case where the user previously had `SPADES_ENABLE_PROJECTS` set to a
#    non-empty list but now the user wishes to disable building all other
#    projects by setting `SPADES_ENABLE_PROJECTS` to an empty string. In that
#    case we still need to set the `SPADES_TOOL_${upper_proj}_BUILD` variables
#    so that we disable building all the projects that were previously enabled.
set(SPADES_ENABLE_PROJECTS_USED OFF CACHE BOOL "")
mark_as_advanced(SPADES_ENABLE_PROJECTS_USED)

if (SPADES_ENABLE_PROJECTS_USED OR NOT SPADES_ENABLE_PROJECTS STREQUAL "")
  set(SPADES_ENABLE_PROJECTS_USED ON CACHE BOOL "" FORCE)
  foreach(proj ${SPADES_KNOWN_PROJECTS} ${SPADES_EXTERNAL_PROJECTS})
    string(TOUPPER "${proj}" upper_proj)
    string(REGEX REPLACE "-" "_" upper_proj ${upper_proj})
    if ("${proj}" IN_LIST SPADES_ENABLE_PROJECTS)
      message(STATUS "${proj} project is enabled")
      set(SHOULD_ENABLE_PROJECT TRUE)
      set(PROJ_DIR ${SPADES_MAIN_PROJ_DIR}/${proj})
      if (NOT EXISTS "${PROJ_DIR}" OR NOT IS_DIRECTORY "${PROJ_DIR}")
        message(FATAL_ERROR "SPADES_ENABLE_PROJECTS requests ${proj} but directory not found: ${PROJ_DIR}")
      endif()
      if (SPADES_EXTERNAL_${upper_proj}_SOURCE_DIR STREQUAL "")
        set(SPADES_EXTERNAL_${upper_proj}_SOURCE_DIR "${PROJ_DIR}" CACHE PATH "" FORCE)
      else()
        set(SPADES_EXTERNAL_${upper_proj}_SOURCE_DIR "${PROJ_DIR}" CACHE PATH "")
      endif()
    elseif ("${proj}" IN_LIST SPADES_EXTERNAL_PROJECTS)
      message(STATUS "${proj} project is enabled")
      set(SHOULD_ENABLE_PROJECT TRUE)
    else()
      message(STATUS "${proj} project is disabled")
      set(SHOULD_ENABLE_PROJECT FALSE)
    endif()
    # Force `SPADES_TOOL_${upper_proj}_BUILD` variables to have values that
    # corresponds with `SPADES_ENABLE_PROJECTS`. This prevents the user setting
    # `SPADES_TOOL_${upper_proj}_BUILD` variables externally.
    set(SPADES_TOOL_${upper_proj}_BUILD
      ${SHOULD_ENABLE_PROJECT}
      CACHE
      BOOL "Whether to build ${upper_proj} as part of SPAdes" FORCE
    )
  endforeach()
endif()
unset(SHOULD_ENABLE_PROJECT)
