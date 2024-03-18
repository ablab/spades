//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#include "kmer_buckets.hpp"
#include "io/binary/binary.hpp"

namespace kmer {

void KMerSegmentPolicyBase::BinRead(std::istream &is) {
    io::binary::BinRead(is, num_segments_);
}

void KMerSegmentPolicyBase::BinWrite(std::ostream &os) const {
    io::binary::BinWrite(os, num_segments_);
}

}
