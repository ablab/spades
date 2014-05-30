############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(ssw C CXX)

add_library(ssw STATIC
            ssw.c
            ssw_cpp.cpp)
            
remove_definitions(-Wall -Wextra -Wconversion -Wno-long-long -Wwrite-strings -pedantic)

target_link_libraries(ssw)

