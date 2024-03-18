############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# This file will recurse into all subdirectories that contain CMakeLists.txt
# Setting variables that match the pattern SPADES_TOOL_{NAME}_BUILD to Off will
# prevent traversing into a directory.
#
# The only tools that need to be explicitly added are ones that have explicit
# ordering requirements.

# Iterates all the subdirectories to create CMake options to enable/disable
# traversing each directory.
create_spades_tool_options()

# Automatically add remaining sub-directories containing a 'CMakeLists.txt'
# file as external projects.
add_spades_implicit_projects()

# Add subprojects specified using SPADES_EXTERNAL_PROJECTS
foreach(p ${SPADES_EXTERNAL_PROJECTS})
  add_spades_external_project(${p})
endforeach(p)
