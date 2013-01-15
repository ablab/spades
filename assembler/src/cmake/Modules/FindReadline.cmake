# GNU Readline library finder
#FIND_PATH(READLINE_INCLUDE_DIR readline/readline.h)
#FIND_LIBRARY(READLINE_LIBRARY NAMES readline) 

#IF (READLINE_INCLUDE_DIR AND READLINE_LIBRARY)
  #SET(READLINE_FOUND TRUE)
  #MESSAGE(STATUS "Found Readline library: ${READLINE_LIBRARY}")
  #MESSAGE(STATUS "Include dir is: ${READLINE_INCLUDE_DIR}")
  #INCLUDE_DIRECTORIES(${READLINE_INCLUDE_DIR})
#ELSE (READLINE_INCLUDE_DIR AND READLINE_LIBRARY)
  #SET(READLINE_FOUND FALSE)
  #MESSAGE(FATAL_ERROR "** Readline library not found!\n** Your distro may provide a binary for Readline e.g. for ubuntu try apt-get install libreadline5-dev")
#ENDIF (READLINE_INCLUDE_DIR AND READLINE_LIBRARY)
include(FindPackageHandleStandardArgs)
if(READLINE_INCLUDE_DIR AND READLINE_LIBRARY AND NCURSES_LIBRARY)
  set(READLINE_FOUND TRUE)
else(READLINE_INCLUDE_DIR AND READLINE_LIBRARY AND NCURSES_LIBRARY)
  FIND_PATH(READLINE_INCLUDE_DIR readline/readline.h /usr/include/readline)

# 2008-04-22 The next clause used to read like this:
#
#  FIND_LIBRARY(READLINE_LIBRARY NAMES readline)
#        FIND_LIBRARY(NCURSES_LIBRARY NAMES ncurses )
#        include(FindPackageHandleStandardArgs)
#        FIND_PACKAGE_HANDLE_STANDARD_ARGS(Readline DEFAULT_MSG NCURSES_LIBRARY READLINE_INCLUDE_DIR READLINE_LIBRARY )
#
# I was advised to modify it such that it will find an ncurses library if
# required, but not if one was explicitly given, that is, it allows the
# default to be overridden. PH 

  FIND_LIBRARY(READLINE_LIBRARY NAMES readline)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(readline DEFAULT_MSG READLINE_INCLUDE_DIR READLINE_LIBRARY )
  mark_as_advanced(READLINE_INCLUDE_DIR READLINE_LIBRARY)
endif(READLINE_INCLUDE_DIR AND READLINE_LIBRARY AND NCURSES_LIBRARY)
