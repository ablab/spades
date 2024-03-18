
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gap_info.hpp"

#include "io/reads/single_read.hpp"

namespace omnigraph {

Sequence Subseq(const io::SingleRead& read, size_t start, size_t end) {
    VERIFY(end > start);
    auto subread = read.Substr(start, end);
    if (subread.IsValid()) {
        return subread.sequence();
    } else {
        return Sequence();
    }
}

Sequence Subseq(const io::SingleReadSeq& read, size_t start, size_t end) {
    return read.sequence().Subseq(start, end);
}

Sequence Subseq(const Sequence& s, size_t start, size_t end) {
    return s.Subseq(start, end);
}

}

