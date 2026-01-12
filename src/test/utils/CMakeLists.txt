############################################################################
# Copyright (c) 2025 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(utils_test CXX)

add_executable(utils_test
               error_codes_test.cpp
               test.cpp)
target_link_libraries(utils_test common_modules ${COMMON_LIBRARIES} teamcity_gtest gtest)

add_test(NAME utils_test COMMAND utils_test)
