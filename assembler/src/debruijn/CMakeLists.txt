project(debruijn CXX)

add_executable(debruijn
 main.cpp
 debruijn_graph.cpp)

target_link_libraries(debruijn ${COMMON_LIBRARIES})

