############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(test_distance CXX)

include_directories(../../../hammer-it)

add_executable(test_distance
              ../../../hammer-it/err_helper_table.cpp
              test_distance.cpp)

target_link_libraries(test_distance ${COMMON_LIBRARIES})
