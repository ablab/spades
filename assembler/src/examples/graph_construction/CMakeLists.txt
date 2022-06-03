set(PROJECT_NAME graph_construction)
PROJECT(${PROJECT_NAME})
add_executable(${PROJECT_NAME} graph_construction.cpp)

target_link_libraries(${PROJECT_NAME}
        utils
        input
        examples_graph_io
        assembly_graph
        mph_index
        sequence
        gqf
        coverage_model)