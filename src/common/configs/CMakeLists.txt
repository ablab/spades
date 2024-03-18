############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(configs CXX)

add_library(configs STATIC config_struct.cpp pe_config_struct.cpp)

target_link_libraries(configs input llvm-support library)