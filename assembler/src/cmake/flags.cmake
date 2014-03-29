# -*- cmake -*-

# Handle OpenMP flags
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# Use parallel libstdc++ if possible
  add_definitions(-DUSE_GLIBCXX_PARALLEL=1)
else ()
  if (NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message(FATAL_ERROR "SPAdes requires OpenMP to be available")
  endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wno-deprecated")
add_definitions(-Wall -Wextra -Wconversion -Wno-sign-conversion -Wno-long-long -Wwrite-strings)
#add_definitions(-Wall)
if (NOT OPENMP_FOUND)
  add_definitions(-Wno-unknown-pragmas)
endif()

# Use libc++ with clang due to C++11 mode
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  # Require libsupc++ on Linux
  if (UNIX AND NOT APPLE)
    set(SYSTEM_LIBRARIES "supc++")
  endif()
endif()

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  message("Making Debug Configuration...")

  add_definitions(-g3)
  add_definitions(-D_GLIBCXX_DEBUG)
  set(SPADES_DEBUG_LOGGING)
else()
  message("Making Release Configuration...")

  if (${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    add_definitions(-g3)
  else()
    add_definitions(-g0)
  endif()

  add_definitions(-O2)
  if (${CMAKE_BUILD_TYPE} STREQUAL "RelWithAsserts" OR
      ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    add_definitions(-UNDEBUG)
  else()
    add_definitions(-DNDEBUG)
  endif()
endif()

# Make sure we're building with frame pointer if tcmalloc is in use
if (SPADES_USE_TCMALLOC)
  add_definitions(-fno-omit-frame-pointer)
endif()

# We will need this at least for jemalloc
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-all_load")
endif()
