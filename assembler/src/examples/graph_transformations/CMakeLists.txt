set(PROJECT_NAME graph_transformations)

project(${PROJECT_NAME})
add_executable(${PROJECT_NAME} graph_transformations.cpp)

target_link_libraries(${PROJECT_NAME}
        examples_graph_io
        utils
        assembly_graph
        #input ?
        )
