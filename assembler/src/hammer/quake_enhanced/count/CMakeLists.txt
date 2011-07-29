project(quake CXX)

add_executable(count
   quake_count.cpp)

target_link_libraries(quake_count ${COMMON_LIBRARIES})




