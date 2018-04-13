############################################################################
# Copyright (c) 2015 Saint Petersburg State University
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
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/corrector"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
