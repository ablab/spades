############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(stages CXX)

add_library(stages STATIC
            construction.cpp simplification.cpp)
target_link_libraries(stages gqf)
