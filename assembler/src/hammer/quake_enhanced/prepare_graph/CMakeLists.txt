project(quake CXX)

add_executable(prepare_graph
   prepare_graph.cpp)

target_link_libraries(prepare_graph ${COMMON_LIBRARIES})




