project(input CXX)

add_library(input
 parser.cpp)

target_link_libraries(input ${COMMON_LIBRARIES})

