cmake_minimum_required(VERSION 2.8)

project(kmer_count_filter CXX)

include_directories(kmc_api)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${EXT_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/debruijn)

add_executable(kmer_multiplicity_counter
        kmc_api/kmc_file.cpp
        kmc_api/kmer_api.cpp
        kmc_api/mmer.cpp
        kmer_multiplicity_counter.cpp)

target_link_libraries(kmer_multiplicity_counter input ${COMMON_LIBRARIES})

add_executable(annotation_propagator
               propagate.cpp)

target_link_libraries(annotation_propagator input mph_index nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})

add_executable(read_binning
               read_binning.cpp)

target_link_libraries(read_binning input mph_index nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})

add_executable(stats
               stats.cpp
               ../debruijn/genome_storage.cpp)
               
target_link_libraries(stats input mph_index nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})            

add_executable(contig_abundance_counter
               contig_abundance_counter.cpp)

target_link_libraries(contig_abundance_counter input mph_index nlopt BamTools ssw yaml-cpp ${COMMON_LIBRARIES})
