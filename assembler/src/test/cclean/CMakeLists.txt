############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cclean_test CXX)

aux_source_directory(. SRC_LIST)
include_directories("${CMAKE_SOURCE_DIR}/cclean")
file(GLOB ${CMAKE_CURRENT_SOURCE_DIR}
    "*.hh"
    "*.h"
    "*.hpp"
    "*.cpp"
)
add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(cclean_test ${COMMON_LIBRARIES} input cityhash otkpp ssw ${ZLIB_LIBRARIES})
