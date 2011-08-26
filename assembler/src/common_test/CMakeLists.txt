project(common_test CXX)

add_executable(common_test
 test.cpp)

target_link_libraries(common_test ${COMMON_LIBRARIES} input)

