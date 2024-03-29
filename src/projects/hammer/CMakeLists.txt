############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(hammer CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(spades-hammer
               main.cpp
               hammer_tools.cpp
               hamcluster.cpp
               kmer_cluster.cpp
               kmer_data.cpp
               config_struct_hammer.cpp
               read_corrector.cpp
               expander.cpp)

target_link_libraries(spades-hammer Blaze common_modules input utils mph_index pipeline gqf ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-hammer PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-hammer
        DESTINATION bin
        COMPONENT spades)
install(DIRECTORY "configs/"
        DESTINATION share/spades/configs/hammer
        COMPONENT spades
        FILES_MATCHING PATTERN "*.info")
