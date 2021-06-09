############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


project(bin-refine CXX)

set(BLAZE_LIBRARY_PATH ${EXT_DIR}/include/blaze)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(bin-refine
               binning_refiner.cpp binning.cpp labels_propagation.cpp link_index.cpp paired_end.cpp
               binning_assignment_strategy.cpp majority_length_strategy.cpp max_likelihood_strategy.cpp)
target_include_directories(bin-refine PRIVATE ${BLAZE_LIBRARY_PATH})
target_link_libraries(bin-refine graphio common_modules ${COMMON_LIBRARIES})

install(TARGETS bin-refine
        DESTINATION bin
        COMPONENT runtime)
