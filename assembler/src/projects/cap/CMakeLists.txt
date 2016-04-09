############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cap CXX)

# Define minimum and maximum K
set(SPADES_MIN_K 1 CACHE INTEGER "Minimum k-mer length")
set(SPADES_MAX_K 300 CACHE INTEGER "Maximum k-mer length")
configure_file("${SPADES_MAIN_INCLUDE_DIR}/k_range.hpp.in"
               "${SPADES_BUILT_INCLUDE_DIR}/k_range.hpp")
include_directories(${CMAKE_SOURCE_DIR}/modules
                    ${CMAKE_SOURCE_DIR}/online_vis)

add_executable(cap-tools
  ${EXT_DIR}/include/teamcity_boost/teamcity_boost.cpp
  ${EXT_DIR}/include/teamcity_boost/teamcity_messages.cpp
  tools.cpp)

add_executable(cap main.cpp)

#if (DEFINED static_build)
#  set_target_properties(cap PROPERTIES LINK_SEARCH_END_STATIC 1)
#endif()

if (READLINE_FOUND)
  include_directories(${READLINE_INCLUDE_DIR})
else(READLINE_FOUND)
  #MESSAGE(WARNING "Library `readline' was not found (not installed?).")
endif()
if (CURSES_FOUND)
  include_directories(${CURSES_INCLUDE_PATH})
else(CURSES_FOUND)
  #MESSAGE(WARNING "Library `ncurses' was not found (not installed?)")
endif()
include_directories(./drawing_commands)

set(CAP_LIBS input cityhash nlopt ${COMMON_LIBRARIES})
if (READLINE_FOUND AND CURSES_FOUND)
  set(CAP_LIBS ${CAP_LIBS} ${READLINE_LIBRARY} ${CURSES_NCURSES_LIBRARY})
endif()

target_link_libraries(cap pipeline dev-support path_extend ${CAP_LIBS})
target_link_libraries(cap-tools dev-support math ${CAP_LIBS})
