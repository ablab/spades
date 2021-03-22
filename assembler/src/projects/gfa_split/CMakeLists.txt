############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


project(gfa-split CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(gfa-split main.cpp)
target_link_libraries(gfa-split toolchain graphio utils assembly_graph ${COMMON_LIBRARIES})

install(TARGETS gfa-split
        DESTINATION bin
        COMPONENT runtime)
