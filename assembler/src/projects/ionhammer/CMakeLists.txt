############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(ionhammer CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(spades-ionhammer
               kmer_data.cpp
               subcluster.cpp
               err_helper_table.cpp
               config_struct.cpp
               kmer_helpers.cpp
               seqeval/BaseHypothesisEvaluator.cpp
               seqeval/TreephaserLite.cpp
               main.cpp
               reference.cpp
               quality_metrics.cpp
               quality_thresholds_estimator.cpp
               hamcluster_1.cpp
               gamma_poisson_model.cpp
               normal_quality_model.cpp)

#add_executable(kmer_evaluator
#               kmer_data.cpp
#               kmer_evaluator.cpp
#               kmer_helpers.cpp
#               err_helper_table.cpp)



target_link_libraries(spades-ionhammer modules input utils pipeline mph_index ${COMMON_LIBRARIES})
#target_link_libraries(kmer_evaluator input  utils mph_index  BamTools ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-ionhammer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-ionhammer
        COMPONENT spades
        RUNTIME DESTINATION bin)
install(DIRECTORY "ionhammer/"
        DESTINATION share/spades/configs/ionhammer
        COMPONENT spades
        FILES_MATCHING PATTERN "*.cfg")
