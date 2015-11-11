cmake_minimum_required(VERSION 2.8)
project(kmer_count_filter)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(kmc_api)
include_directories(../include)
include_directories(../../ext/include)

add_executable(kmer_multiplicity_counter
        kmc_api/kmc_file.cpp
        kmc_api/kmer_api.cpp
        kmc_api/mmer.cpp
        kmer_multiplicity_counter.cpp)
