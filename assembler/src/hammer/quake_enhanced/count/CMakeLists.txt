############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(quake CXX)

add_executable(count
   count.cpp)

target_link_libraries(count ${COMMON_LIBRARIES})




