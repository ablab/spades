############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(gzstream CXX)

add_library(gzstream STATIC gzstream.C)
include_directories(SYSTEM "${ZLIB_INCLUDE_DIRS}")
target_link_libraries(gzstream ${ZLIB_LIBRARIES})
