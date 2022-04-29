set(PROJECT_NAME examples_graph_io)

project(${PROJECT_NAME})
add_library(${PROJECT_NAME} gfa_io.cpp)

target_link_libraries(${PROJECT_NAME} graphio assembly_graph)