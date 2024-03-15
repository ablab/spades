############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(include_test CXX)

add_executable(include_test
               seq_test.cpp sequence_test.cpp rtseq_test.cpp quality_test.cpp nucl_test.cpp
               cyclic_hash_test.cpp binary_test.cpp
               test.cpp)
target_link_libraries(include_test common_modules input ${COMMON_LIBRARIES} teamcity_gtest gtest)

add_test(NAME include_test COMMAND include_test)

