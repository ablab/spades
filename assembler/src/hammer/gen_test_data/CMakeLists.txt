############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(gen_test_data CXX)

add_executable(gen_test_data
   main.cpp)

target_link_libraries(gen_test_data ${COMMON_LIBRARIES})




