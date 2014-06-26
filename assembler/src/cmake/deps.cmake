# -*- cmake -*-

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  # Require at least gcc 4.7
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
    message(FATAL_ERROR "SPAdes requires gcc version 4.7 or later")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
    message(FATAL_ERROR "SPAdes requires clang version 3.2 or later")
  endif()
else()
  message(WARNING "Unsupported compiler is detected. SPAdes compilation was not tested on it and may fail")
endif()

find_package(OpenMP)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(Readline QUIET)
set(CURSES_NEED_NCURSES TRUE)
find_package(Curses QUIET)

# Use included boost unless explicitly specified
if (NOT SPADES_BOOST_ROOT)
  set(BOOST_ROOT "${EXT_DIR}/include")
else()
  set(BOOST_ROOT SPADES_BOOST_ROOT)
endif()
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED)
