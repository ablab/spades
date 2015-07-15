############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(online_vis CXX)

add_executable(online_vis
               ${CMAKE_SOURCE_DIR}/debruijn/kmer_coverage_model.cpp
               ${CMAKE_SOURCE_DIR}/debruijn/config_struct.cpp
               ${CMAKE_SOURCE_DIR}/debruijn/stage.cpp
               ${CMAKE_SOURCE_DIR}/debruijn/genomic_info_filler.cpp
               ${CMAKE_SOURCE_DIR}/debruijn/path_extend/pe_config_struct.cpp
               autocompletion.cpp
               main.cpp)

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
include_directories(${CMAKE_SOURCE_DIR}/debruijn)

if (READLINE_FOUND AND CURSES_FOUND)
  target_link_libraries(online_vis input nlopt mph_index format ${COMMON_LIBRARIES} ${READLINE_LIBRARY} ${CURSES_NCURSES_LIBRARY})
elseif (READLINE_FOUND)
  target_link_libraries(online_vis input nlopt mph_index format ${COMMON_LIBRARIES} ${READLINE_LIBRARY})
else()
  target_link_libraries(online_vis input nlopt mph_index format ${COMMON_LIBRARIES})
endif()

if (DEFINED static_build)
  set_target_properties(online_vis PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

