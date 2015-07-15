############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(quake CXX)

add_executable(filter_trusted_enh
   main.cpp)

target_link_libraries(filter_trusted_enh ${COMMON_LIBRARIES})
