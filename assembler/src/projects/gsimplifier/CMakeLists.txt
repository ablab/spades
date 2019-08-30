############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# Copyright (c) 2019 University of Warwick
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades-gsimplifier CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(spades-gsimplifier
               ../unitig_coverage/profile_storage.cpp
               main.cpp)

target_link_libraries(spades-gsimplifier graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-gsimplifier PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-gsimplifier
        DESTINATION bin
        COMPONENT runtime)
