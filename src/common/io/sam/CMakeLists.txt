############################################################################
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(samio CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(samio STATIC
            read.cpp sam_reader.cpp bam_parser.cpp)
target_link_libraries(samio BamTools samtools ${ZLIB_LIBRARIES})
