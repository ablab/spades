############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2020-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(toolchain CXX)

add_library(toolchain STATIC
            utils.cpp)
target_link_libraries(toolchain graphio binary_io utils)
