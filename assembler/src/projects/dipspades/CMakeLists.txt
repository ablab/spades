############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(dipspades CXX)

add_executable(dipspades
               dipspades_config.cpp
               utils/files_utils.cpp
               main.cpp)

target_link_libraries(dipspades common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(dipspades PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS dipspades
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/dipspades"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
