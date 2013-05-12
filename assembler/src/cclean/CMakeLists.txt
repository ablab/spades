############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cclean CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(cclean
               aho_corasick.cpp
               database.cpp
               utils.cpp
               QcException.cpp
               running_modes.cpp
               job_wrappers.cpp
               output.cpp
               config_struct_cclean.cpp
               adapter_index.cpp
               main.cpp)

target_link_libraries(cclean ssw input mph_index ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(cclean PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS cclean
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/cclean"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
