############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(quake_count CXX)

add_executable(quake_count
   quake_count.cpp)

target_link_libraries(quake_count ${COMMON_LIBRARIES})




