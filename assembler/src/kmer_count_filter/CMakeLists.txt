cmake_minimum_required(VERSION 2.8)
project(kmer_count_filter)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(kmc_api)
include_directories(../include)
include_directories(../../ext/include)

add_executable(kmer_count_filter
        kmc_api/kmc_file.cpp
        kmc_api/kmer_api.cpp
        kmc_api/mmer.cpp
        main.cpp)
