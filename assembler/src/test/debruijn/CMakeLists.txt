project(debruijn_test CXX)

add_executable(debruijn_test
 test.cpp)

target_link_libraries(debruijn_test ${COMMON_LIBRARIES} input)

