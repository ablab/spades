############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(long_contigs CXX)

add_executable(long_contigs
 main.cpp)

target_link_libraries(long_contigs ${COMMON_LIBRARIES} input)

