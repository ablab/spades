# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# src/api/
# ==========================

# list include paths
include_directories( ${samtools_SOURCE_DIR}/src )

# add compiler definitions 
add_definitions( -DBAMTOOLS_API_LIBRARY ) # (for proper exporting of library symbols)

# fetch all internal source files

# make list of all API source files
file(GLOB sources "[a-zA-Z]*.c" "bcftools/[a-zA-Z]*.c" "misc/[a-z][A-Z]*.c")



