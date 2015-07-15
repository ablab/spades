############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(quake CXX)

add_executable(test_correction_quality
   main.cpp)

target_link_libraries(test_correction_quality ${COMMON_LIBRARIES})
