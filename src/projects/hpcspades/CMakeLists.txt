############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(hpcspades CXX)

add_subdirectory(common)

add_library(spades-stages-hpc STATIC distance_estimation_mpi.cpp)

set(HPCSPADES_SRC pipeline.cpp ../../projects/spades/series_analysis.cpp ../../projects/mts/contig_abundance.cpp)
set(HPCSPADES_LIB spades-stages-hpc spades-stages graphio common_modules paired_info_mpi stages-mpi ${COMMON_LIBRARIES})

add_executable(spades-hpc main_mpi.cpp ${HPCSPADES_SRC})

target_link_libraries(spades-hpc ${HPCSPADES_LIB} ${MPI_LIBRARIES})
set_target_properties(spades-hpc PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

if (SPADES_STATIC_BUILD)
    set_target_properties(spades-hpc PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS spades-hpc
        DESTINATION bin
        COMPONENT runtime)
