############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(quake_correct CXX)

add_executable(correct
   bithash.cpp correct.cpp edit.cpp Read.cpp gzstream.C)

target_link_libraries(correct ${COMMON_LIBRARIES})




