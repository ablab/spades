############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


project(bin-refine CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(bin-refine binning_refiner.cpp binning.cpp)
target_link_libraries(bin-refine graphio common_modules ${COMMON_LIBRARIES})

install(TARGETS bin-refine
        DESTINATION bin
        COMPONENT runtime)
