project(input CXX)

add_library(input STATIC
 parser.cpp)

target_link_libraries(input ${COMMON_LIBRARIES})

