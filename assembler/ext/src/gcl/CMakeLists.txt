############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(gcl CXX)

include_directories(${EXT_DIR}/include/gcl)

add_library(gcl STATIC
    atomic.c)
