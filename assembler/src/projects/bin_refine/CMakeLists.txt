############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


project(bin-refine CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(bin-refine
               alpha_assigner.cpp alpha_propagation.cpp
               binning_refiner.cpp binning.cpp labels_propagation.cpp link_index.cpp
               paired_end.cpp read_splitting.cpp
               binning_assignment_strategy.cpp majority_length_strategy.cpp max_likelihood_strategy.cpp)
target_link_libraries(bin-refine Blaze graphio toolchain common_modules ${COMMON_LIBRARIES})

install(TARGETS bin-refine
        DESTINATION bin
        COMPONENT runtime)
