############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(stages CXX)

set(stages_src read_conversion.cpp construction.cpp simplification.cpp ss_edge_split.cpp genomic_info_filler.cpp)
add_library(stages STATIC ${stages_src})

target_link_libraries(stages coverage_model pipeline gqf input)
