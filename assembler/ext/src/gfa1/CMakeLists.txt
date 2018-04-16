############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(gfa1 C)


# Create libraries.
add_library(gfa1 STATIC gfa.c)
target_link_libraries(gfa1 ${ZLIB_LIBRARIES})


