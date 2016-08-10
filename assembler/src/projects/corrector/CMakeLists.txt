############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(corrector CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(corrector
	    positional_read.cpp
        interesting_pos_processor.cpp
        contig_processor.cpp
        dataset_processor.cpp
        config_struct.cpp
        main.cpp)

target_link_libraries(corrector input common_modules ${COMMON_LIBRARIES})



if (SPADES_STATIC_BUILD)
  set_target_properties(corrector PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS corrector
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/corrector"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
