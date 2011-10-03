project(utils CXX)

add_executable(utils
 paired_info_checker.cpp)

target_link_libraries(utils ${COMMON_LIBRARIES} input)