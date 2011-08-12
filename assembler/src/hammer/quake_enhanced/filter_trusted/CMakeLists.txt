project(quake CXX)

add_executable(filter_trusted
   main.cpp)

target_link_libraries(filter_trusted ${COMMON_LIBRARIES})
