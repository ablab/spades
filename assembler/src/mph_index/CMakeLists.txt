############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(mph_index CXX)

add_library(mph_index STATIC
            MurmurHash3.cpp
            mph_bits.cc
            mph_index.cc
            trigraph.cc)
