//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <io/sam/read.hpp>

using namespace std;

namespace sam_reader {

string SingleSamRead::cigar() const {
    uint32_t *cigar = bam1_cigar(data_);
    string res;
    res.reserve(data_->core.n_cigar);
    for (size_t k = 0; k < data_->core.n_cigar; ++k) {
        res += std::to_string(bam_cigar_oplen(cigar[k]));
        res += bam_cigar_opchr(cigar[k]);

    }
    return res;
}

string SingleSamRead::name() const {
    string res(bam1_qname(data_));
    return res;
}

string SingleSamRead::seq() const {
    string res = "";
    auto b = bam1_seq(data_);
    for (int k = 0; k < data_->core.l_qseq; ++k) {
        res += bam_nt16_rev_table[bam1_seqi(b, k)];
    }
    return res;
}


}
;
