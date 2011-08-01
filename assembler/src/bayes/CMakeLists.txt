project(bayes CXX)

find_package(OpenMP)

# bug in cmake, fixed in 2.6.4 (patch 8)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")

add_definitions(${OpenMP_CXX_FLAGS})

add_executable(bayes
   main.cpp
   bayes_quality.cpp)

target_link_libraries(bayes ${COMMON_LIBRARIES})

