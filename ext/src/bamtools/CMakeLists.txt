# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# top-level
# ==========================

# set project name
project( BamTools )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include)

# set BamTools version information
set( BamTools_VERSION_MAJOR 2 )
set( BamTools_VERSION_MINOR 3 )
set( BamTools_VERSION_BUILD 0 )

add_definitions( -D_FILE_OFFSET_BITS=64 )

# list subdirectories to build in
add_subdirectory( api )
