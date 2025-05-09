############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(utils CXX)

set(utils_src
    memory_limit.cpp
    filesystem/path_helper.cpp
    filesystem/temporary.cpp
    filesystem/glob.cpp
    logger/logger_impl.cpp
    logger/log_writers.cpp
    logger/log_writers_thread.cpp
    )

if (READLINE_FOUND)
  set(utils_src ${utils_src} autocompletion.cpp)
endif()

add_library(utils STATIC
            ${utils_src})
target_link_libraries(utils ${COMMON_LIBRARIES})
# This is hack, but otherwise it is very hard to obtain additional library paths exposed via FindOpenMP
if (OPENMP_FOUND)
  target_link_libraries(utils OpenMP::OpenMP_CXX)
endif()

add_library(version STATIC
            version.cpp)
