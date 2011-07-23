project(hammer CXX)

find_package(OpenMP)

# bug in cmake, fixed in 2.6.4 (patch 8)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")

add_definitions(${OpenMP_CXX_FLAGS})
add_definitions(-fPIC)

add_executable(main
   main.cpp
   hammer_tools.cpp
   kmer_cluster.cpp)
   
add_executable(prepare_graph
   prepare_graph.cpp)
   
add_executable(preproc
   preproc.cpp)

target_link_libraries(main ${COMMON_LIBRARIES})
target_link_libraries(prepare_graph ${COMMON_LIBRARIES})
target_link_libraries(preproc ${COMMON_LIBRARIES})


