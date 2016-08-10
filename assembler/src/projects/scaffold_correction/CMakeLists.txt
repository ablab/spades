############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(moleculo CXX)


add_executable(scaffold_correction
               main.cpp)
target_link_libraries(scaffold_correction common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(scaffold_correction PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS scaffold_correction
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/scaffold_correction"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
