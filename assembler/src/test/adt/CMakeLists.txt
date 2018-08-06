############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(adt_test CXX)

add_executable(pack_demo
               pack_demo.cpp)
target_link_libraries(pack_demo ${COMMON_LIBRARIES})
