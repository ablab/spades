project(kmer_count CXX)

add_executable(kmer_count
 kmer_count.cpp)

target_link_libraries(kmer_count ${COMMON_LIBRARIES} input)