############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(utils CXX)

add_executable(utils
 paired_info_checker.cpp)

target_link_libraries(utils ${COMMON_LIBRARIES} input)