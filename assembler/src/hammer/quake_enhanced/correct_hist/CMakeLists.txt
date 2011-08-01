project(quake CXX)

add_executable(correct_hist
   main.cpp)

target_link_libraries(correct_hist ${COMMON_LIBRARIES})
