//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * read.hpp
 *
 *  Created on: Jun 26, 2014
 *      Author: lab42
 */
#include "samtools/bam.h"

#include <string>
#include <unordered_map>
#include <samtools/bam.h>

#pragma once

namespace sam_reader {

class SingleSamRead {
private:
    bam1_t *data_;

public:

    SingleSamRead() {
        data_ = bam_init1();
    }
    SingleSamRead(SingleSamRead const &c) {
        data_ = bam_dup1( c.data_);
    }
    ~SingleSamRead() {
        bam_destroy1(data_);
    }
    SingleSamRead& operator= (const SingleSamRead &c){
        bam_destroy1(data_);
        data_ = bam_dup1(c.data_);
        return *this;
    }

    int32_t get_data_len() const {
        return data_->core.l_qseq;
    }

    uint32_t get_cigar_len() const {
        return data_->core.n_cigar;
    }

    int get_contig_id() const {
        return data_->core.tid;
    }

    bool is_aligned() const {
        return (data_->core.flag & 0x4) == 0;
    }

    bool is_properly_aligned() const {
        return is_aligned() && get_map_qual() != 0;
    }

    bool get_strand() const {
        return (data_->core.flag & 0x10) == 0;
    }

    uint32_t get_map_qual() const {
        return data_->core.qual;
    }

    int32_t get_pos() const {
        return data_->core.pos;
    }

    uint32_t* get_cigar_ptr() const {
        return bam1_cigar(data_);
    }

    uint8_t* get_seq_ptr() const {
        return bam1_seq(data_);
    }

    std::string get_cigar() const;
    std::string get_name() const;
    std::string get_seq() const;

    void set_data(bam1_t *seq_) {
        bam_destroy1(data_);
        data_ = bam_dup1( seq_);
    }
};

class PairedSamRead {
private:
    SingleSamRead r1;
    SingleSamRead r2;

public:
    PairedSamRead(): r1(), r2() {
    }

    PairedSamRead(SingleSamRead &a1, SingleSamRead &a2) {
        r1 = a1;
        r2 = a2;
    }

    const SingleSamRead& GetLeft() const {
        return r1;
    }

    const SingleSamRead& GetRight() const {
        return r2;
    }
};
}
;
