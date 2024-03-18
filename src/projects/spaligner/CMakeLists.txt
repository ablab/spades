############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spaligner CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(spaligner
               align_longreads.cpp
               mapping_printer.cpp)
target_link_libraries(spaligner common_modules bwa edlib graphio ${COMMON_LIBRARIES})

add_executable(form_truealignments form_truealignments.cpp)
target_link_libraries(form_truealignments common_modules bwa edlib graphio ${COMMON_LIBRARIES})

install(TARGETS spaligner
        DESTINATION bin
        COMPONENT spaligner)

install(FILES spaligner_config.yaml
        COMPONENT spaligner
        DESTINATION share/spaligner)
