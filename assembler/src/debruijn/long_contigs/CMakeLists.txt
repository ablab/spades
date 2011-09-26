project(long_contigs CXX)

add_executable(long_contigs
 main.cpp)

target_link_libraries(long_contigs ${COMMON_LIBRARIES} input)

