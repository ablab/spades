############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2021-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


project(binspreader CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(binspreader
               alpha_assigner.cpp alpha_propagation.cpp
               binning_refiner.cpp binning.cpp labels_propagation.cpp link_index.cpp
               paired_end.cpp read_splitting.cpp
               binning_assignment_strategy.cpp majority_length_strategy.cpp max_likelihood_strategy.cpp)
target_link_libraries(binspreader Blaze graphio toolchain common_modules ${COMMON_LIBRARIES})

install(TARGETS binspreader
        DESTINATION bin
        COMPONENT binspreader)
