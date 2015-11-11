cmake_minimum_required(VERSION 2.8)

project(kmer_count_filter CXX)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(kmc_api)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${EXT_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/debruijn)

add_executable(kmer_multiplicity_counter
        kmc_api/kmc_file.cpp
        kmc_api/kmer_api.cpp
        kmc_api/mmer.cpp
        kmer_multiplicity_counter.cpp)

add_executable(annotation_propagator
               propagate.cpp)

target_link_libraries(annotation_propagator input mph_index nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})

