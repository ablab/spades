/*
 * read.hpp
 *
 *  Created on: Jun 26, 2014
 *      Author: lab42
 */
#include "positional_read.hpp"

#include "samtools/bam.h"

#include <string>
#include <unordered_map>

#pragma once

namespace corrector {

struct SingleSamRead {
    bam1_t *data_;
    size_t get_data_len() const {
        return data_->core.l_qseq;
    }
    size_t get_cigar_len() const {
        return data_->core.n_cigar;
    }
    int get_contig_id() const {
        return data_->core.tid;
    }
    void set_data(bam1_t *seq_) {
        bam_destroy1(data_);
        data_ = bam_dup1( seq_);
    }

    std::string get_cigar() const;
    std::string get_name() const;
    std::string get_seq() const;

    int CountPositions(std::unordered_map<size_t, position_description> &ps, const size_t &contig_length) const;
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
};

struct PairedSamRead {
    SingleSamRead r1;
    SingleSamRead r2;
    PairedSamRead(): r1(), r2(){
    }

    PairedSamRead(SingleSamRead &a1, SingleSamRead &a2) {
        r1 = a1;
        r2 = a2;
    }
    int CountPositions(std::unordered_map<size_t, position_description> &ps, const size_t &contig_length) const;
};
}
;
