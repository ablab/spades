############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(utils CXX)

set(utils_src
    copy_file.cpp
    path_helper.cpp
    logger/logger_impl.cpp)

if (READLINE_FOUND)
  set(utils_src ${utils_src} autocompletion.cpp)
endif()

add_library(utils STATIC
            ${utils_src})
