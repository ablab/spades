project(iolib CXX)

add_library(iolib
 parser.hpp
 parser.cpp)

target_link_libraries(iolib ${COMMON_LIBRARIES})

