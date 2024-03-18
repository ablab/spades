############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2018-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(mph_index CXX)

add_library(mph_index STATIC kmer_mph/kmer_buckets.cpp)
