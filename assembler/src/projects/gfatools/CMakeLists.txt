set(PROJECT_NAME gfa_tools)
project(${PROJECT_NAME})

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} assembly_graph graphio utils)
