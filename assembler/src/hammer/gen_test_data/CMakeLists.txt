project(gen_test_data CXX)

add_executable(gen_test_data
   main.cpp)

target_link_libraries(gen_test_data ${COMMON_LIBRARIES})




