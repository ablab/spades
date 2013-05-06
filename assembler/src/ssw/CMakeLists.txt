############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(ssw C CXX)

add_library(ssw STATIC
            ssw.c
            ssw_cpp.cpp)

target_link_libraries(ssw)

