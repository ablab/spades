############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(galigner CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_executable(longreads_aligner align_longreads.cpp)
target_link_libraries(longreads_aligner common_modules cityhash bwa edlib graphio ${COMMON_LIBRARIES})