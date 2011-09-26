project(utils CXX)

add_executable(paired_info_checker
 paired_info_checker.cpp)

target_link_libraries(paired_info_checker ${COMMON_LIBRARIES} input)