############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(dev_support CXX)

add_library(dev_support STATIC
            copy_file.cpp
            path_helper.cpp
            logger/logger_impl.cpp)

