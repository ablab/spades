# -*- cmake -*-

# Default configuration
set(SPADES_DEFAULT_BUILD_TYPE "RelWithDebInfo" CACHE STRING "SPAdes default build type")
if (NOT CMAKE_BUILD_TYPE)
  message("Setting default build configuration: ${SPADES_DEFAULT_BUILD_TYPE}")
  set(CMAKE_BUILD_TYPE "${SPADES_DEFAULT_BUILD_TYPE}" CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithAsserts RelWithDebInfo."
      FORCE)
endif()

# Define option for turning on/off debug logging
option(SPADES_DEBUG_LOGGING "Turn on debug / trace logging" ON)

option(SPADES_ENABLE_EXPENSIVE_CHECKS "Turn on expensive checks in hot places" OFF)

# Define option to enable / disable ASAN
option(SPADES_ENABLE_ASAN "Turn on / off address sanitizer" OFF)
if (SPADES_ENABLE_ASAN)
  add_compile_options(-fsanitize=address)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

# Define option for static / dynamic build.
option(SPADES_STATIC_BUILD "Link SPAdes statically" OFF)
if (SPADES_STATIC_BUILD)
  # it'll make cmake to find libraries archives, not dynamic link
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a) 
  set(LINK_SEARCH_START_STATIC TRUE)
  set(LINK_SEARCH_END_STATIC TRUE)
  # This is dirty hack to get rid of -Wl,-Bdynamic
  set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")

  if (APPLE)
    message(WARNING "Static build is not quite supported on Mac OS")
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
    endif()
  else()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    add_compile_options(-static)
  endif()

  set(Boost_USE_STATIC_LIBS        ON)
  set(Boost_USE_STATIC_RUNTIME     ON)
endif()

option(SPADES_USE_GPROF "gprof profiler" OFF)

if (SPADES_USE_GPROF)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg") 
endif()

# Define minimum and maximum K
set(SPADES_MIN_K "1" CACHE STRING "Minimum k-mer length")
set(SPADES_MAX_K "128" CACHE STRING "Maximum k-mer length")
configure_file("${SPADES_MAIN_INCLUDE_DIR}/k_range.hpp.in"
               "${SPADES_BUILT_INCLUDE_DIR}/k_range.hpp")

# Define boost root
set(SPADES_BOOST_ROOT "" CACHE PATH "Boost root used to build SPAdes")

# Various internal stuff
option(SPADES_BUILD_INTERNAL "Build internal projects" OFF)
option(SPADES_USE_TCMALLOC "Link spades with TCMalloc" OFF)
if (SPADES_USE_TCMALLOC)
  find_package(GooglePerfTools REQUIRED)

  if (GOOGLE_PERFTOOLS_ROOT)
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif()
endif()
option(SPADES_USE_JEMALLOC "Link spades with JEMalloc" OFF)
if (SPADES_USE_JEMALLOC)
  if (POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
  endif()
  if (SPADES_USE_TCMALLOC OR SPADES_USE_MIMALLOC)
    message(FATAL_ERROR "Cannot link to both TCMalloc/mimalloc and JEMalloc")
  endif()
  set(JEMALLOC_XMALLOC 1 CACHE INTERNAL "" FORCE)
  set(JEMALLOC_BUILD_STATIC ON CACHE INTERNAL "" FORCE)
  set(JEMALLOC_BUILD_SHARED OFF CACHE INTERNAL "" FORCE)
endif()
option(SPADES_USE_MIMALLOC "Link spades with mimalloc" ON)
if (SPADES_USE_MIMALLOC)
  if (SPADES_USE_TCMALLOC OR SPADES_USE_JEMALLOC)
    message(FATAL_ERROR "Cannot link to both TCMalloc/JEMalloc and mimalloc")
  endif()

  set(MI_OSX_ZONE     ON  CACHE INTERNAL "" FORCE)
  set(MI_XMALLOC      ON  CACHE INTERNAL "" FORCE)
  set(MI_BUILD_TESTS  OFF CACHE INTERNAL "" FORCE)
  set(MI_SHOW_ERRORS  ON  CACHE INTERNAL "" FORCE)
  set(MI_BUILD_STATIC OFF CACHE INTERNAL "" FORCE)
  set(MI_BUILD_SHARED OFF CACHE INTERNAL "" FORCE)
  set(MI_BUILD_OBJECT ON  CACHE INTERNAL "" FORCE)
  set(MI_USE_CXX      ON  CACHE INTERNAL "" FORCE)
  if (SPADES_STATIC_BUILD AND NOT APPLE)
    set(MI_OVERRIDE_WRAP ON CACHE INTERNAL "" FORCE)
    set(MI_OVERRIDE_GLIBC OFF CACHE INTERNAL "" FORCE)
  endif()
endif()

option(SPADES_FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." OFF)
if (SPADES_FORCE_COLORED_OUTPUT)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
       add_compile_options(-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       add_compile_options(-fcolor-diagnostics)
    endif()
endif()
