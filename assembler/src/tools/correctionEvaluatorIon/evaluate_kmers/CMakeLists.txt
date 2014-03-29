project(kmer_evaluator CXX)

include_directories(../../../hammer-it)

add_executable(kmer_evaluator
               ../../../hammer-it/kmer_data.cpp
               ../../../hammer-it/err_helper_table.cpp
               kmer_evaluator.cpp)

target_link_libraries(kmer_evaluator input mph_index yaml-cpp ${COMMON_LIBRARIES})
