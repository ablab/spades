############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2018-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(hmmercpp CXX)

add_library(hmmercpp STATIC
            hmmfile.cpp hmmmatcher.cpp)
target_link_libraries(hmmercpp easel hmmer)
