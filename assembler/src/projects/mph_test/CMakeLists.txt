############################################################################
# Copyright (c) 2016 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades-kmercount CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(spades-kmercount
               main.cpp)

target_link_libraries(spades-kmercount spades_modules ${COMMON_LIBRARIES})
