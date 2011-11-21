project(debruijn_test CXX)

include_directories("${CMAKE_SOURCE_DIR}/debruijn")

add_executable(debruijn_test
 test.cpp)

target_link_libraries(debruijn_test ${COMMON_LIBRARIES} input)

