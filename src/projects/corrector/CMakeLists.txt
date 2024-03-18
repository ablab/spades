############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(corrector CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(spades-corrector-core
	      positional_read.cpp
              interesting_pos_processor.cpp
              contig_processor.cpp
              dataset_processor.cpp
              config_struct.cpp
              main.cpp)
target_link_libraries(spades-corrector-core input common_modules ${COMMON_LIBRARIES})



if (SPADES_STATIC_BUILD)
  set_target_properties(spades-corrector-core PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-corrector-core
        DESTINATION bin
        COMPONENT spades)
install(DIRECTORY "configs/"
        DESTINATION share/spades/configs/corrector
        COMPONENT spades
        FILES_MATCHING PATTERN "*.info")
