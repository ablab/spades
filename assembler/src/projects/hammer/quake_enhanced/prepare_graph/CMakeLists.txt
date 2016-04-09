############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(quake CXX)

add_executable(prepare_graph
   prepare_graph.cpp)

target_link_libraries(prepare_graph ${COMMON_LIBRARIES})




