# -*- cmake -*-

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  # Require at least gcc 9.1
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    message(FATAL_ERROR "SPAdes requires gcc version 9.1 or later")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
    message(FATAL_ERROR "SPAdes requires clang version 5 or later")
  endif()
else()
  message(WARNING "Unsupported compiler is detected. SPAdes compilation was not tested on it and may fail")
endif()

find_package(OpenMP COMPONENTS CXX)
find_package(BZip2 REQUIRED)
find_package(Readline QUIET)
set(CURSES_NEED_NCURSES TRUE)
find_package(Curses QUIET)

# See if we can find zstd via cmake config
find_package(zstd QUIET)
if (zstd_FOUND)
  set(SPADES_USE_ZSTD ON)
  if (TARGET zstd::libzstd_static)
    set(ZSTD_LIB zstd::libzstd_static)
  else()
    set(ZSTD_LIB zstd::libzstd_shared)
  endif()
else()
  # No luck, let's try pkg
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(ZSTD REQUIRED libzstd IMPORTED_TARGET)
  if (LIBZSTD_FOUND)
    set(SPADES_USE_ZSTD ON)
    set(ZSTD_LIB PkgConfig::ZSTD)
  endif()
endif()

# Use included boost unless explicitly specified
if (NOT SPADES_BOOST_ROOT)
  set(BOOST_ROOT "${EXT_DIR}/include")
else()
  set(BOOST_ROOT SPADES_BOOST_ROOT)
endif()
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED)
