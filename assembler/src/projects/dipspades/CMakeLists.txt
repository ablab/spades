############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(dipspades CXX)

add_executable(spades-dipspades-core
               dipspades_config.cpp
               utils/files_utils.cpp
               main.cpp)

target_link_libraries(spades-dipspades-core common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-dipspades-core PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-dipspades-core
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/dipspades"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
