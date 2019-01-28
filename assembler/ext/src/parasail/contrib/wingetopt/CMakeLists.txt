PROJECT(wingetopt)
cmake_minimum_required(VERSION 2.8)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_library(wingetopt STATIC src/getopt.c src/getopt.h)

install(FILES src/getopt.h DESTINATION include)
install(TARGETS wingetopt LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

