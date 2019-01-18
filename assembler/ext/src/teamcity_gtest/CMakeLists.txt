############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(teamcity_gtest CXX)

add_library(teamcity_gtest STATIC
            teamcity_gtest.cpp
            teamcity_messages.cpp)
target_link_libraries(teamcity_gtest gtest)

