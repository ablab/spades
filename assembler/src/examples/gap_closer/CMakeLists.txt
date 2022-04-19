set(PROJECT_NAME gap_closer_example)
project(${PROJECT_NAME})

add_executable(${PROJECT_NAME} gap_closer_example.cpp)

target_link_libraries(${PROJECT_NAME}
                        utils
                        common_modules
                        pipeline
                        spades-stages
                        assembly_graph
                        graphio
                        teamcity_gtest gtest)
