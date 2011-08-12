add_subdirectory(count)
add_subdirectory(prepare_graph)
add_subdirectory(filter_trusted)
add_subdirectory(filter_trusted_enh)
add_subdirectory(test_correction_quality)
add_subdirectory(correct_hist)
add_subdirectory(generate_limits)
project(quake_enhanced CXX)
add_executable(quake_enhanced
  main.cpp options.cpp)
target_link_libraries(quake_enhanced ${COMMON_LIBRARIES})



