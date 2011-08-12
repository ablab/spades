project(quake CXX)

add_executable(test_correction_quality
   main.cpp)

target_link_libraries(test_correction_quality ${COMMON_LIBRARIES})
