############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(library CXX)

add_library(library STATIC
            library.cpp
            library_data.cpp)

target_link_libraries(library input llvm-support)

