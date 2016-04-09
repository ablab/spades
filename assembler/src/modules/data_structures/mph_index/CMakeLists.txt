############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(mph_index CXX)

add_library(mph_index STATIC bitpair_vector.cpp)

target_link_libraries(mph_index cityhash)

