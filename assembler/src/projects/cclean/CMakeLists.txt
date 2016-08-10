############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cclean CXX)
aux_source_directory(. SRC_LIST)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB ${CMAKE_CURRENT_SOURCE_DIR}
    "*.hh"
    "*.h"
    "*.hpp"
    "*.cpp"
)
add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(cclean ssw input cityhash ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(cclean PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS cclean
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/cclean"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info")
