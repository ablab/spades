project(quake CXX)

add_executable(generate_limits
   main.cpp)

target_link_libraries(generate_limits ${COMMON_LIBRARIES})
