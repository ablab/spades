############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(pipeline CXX)

add_library(pipeline STATIC
            stage.cpp config_struct.cpp genomic_info_filler.cpp library.cpp)

target_link_libraries(pipeline input llvm-support)

