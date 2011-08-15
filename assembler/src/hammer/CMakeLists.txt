project(hammer CXX C)

find_package(OpenMP)

# bug in cmake, fixed in 2.6.4 (patch 8)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")

add_definitions(${OpenMP_CXX_FLAGS})
add_definitions(-fPIC)

add_executable(main
   main.cpp
   hammer_tools.cpp
   kmer_cluster.cpp
   position_read.cpp
   position_kmer.cpp
   subkmers.cpp)
   
add_subdirectory(quake_count)
add_subdirectory(quake_enhanced)

target_link_libraries(main ${COMMON_LIBRARIES})




