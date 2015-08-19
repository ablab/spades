project(test_distance CXX)

include_directories(../../../hammer-it)

add_executable(test_distance
              ../../../hammer-it/err_helper_table.cpp
              test_distance.cpp)

target_link_libraries(test_distance ${COMMON_LIBRARIES})
