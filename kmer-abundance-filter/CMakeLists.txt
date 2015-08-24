cmake_minimum_required(VERSION 3.2 )
project(kmer_abundance_filter)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(kmer_abundance_filter ${SOURCE_FILES})