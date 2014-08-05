/*
 * read.hpp
 *
 *  Created on: Jun 26, 2014
 *      Author: lab42
 */
#include "include.hpp"
#include "positional_read.hpp"

#include "samtools/bam.h"

#pragma once

namespace corrector {

// WTF: Make sure getters and setters are properly named
struct SingleSamRead {
    bam1_t data_;
    size_t get_data_len() const {
        return data_.core.l_qseq;
    }
    size_t get_cigar_len() const {
        return data_.core.n_cigar;
    }
    int get_contig_id() const {
        return data_.core.tid;
    }
    void set_data(bam1_t *seq_) {
        //TODO: delete
        // WTF: fix TODO
        bam1_t *new_seq = bam_dup1(seq_);
        data_ = *new_seq;
    }
    // WTF: This does not belong here
    int CountPositions(unordered_map<size_t, position_description> &ps, const string &contig) const;
    string get_cigar() const;
    string get_qual() const;
    string get_name() const;
    string get_seq() const;
    ~SingleSamRead() {
    }
};

struct PairedSamRead {
    SingleSamRead r1;
    SingleSamRead r2;
    void pair(SingleSamRead &a1, SingleSamRead &a2);
    PairedSamRead() {
    }
    PairedSamRead(SingleSamRead &a1, SingleSamRead &a2)
            : r1(a1),
              r2(a2) {

    }
    // WTF: This does not belong here    
    int CountPositions(unordered_map<size_t, position_description> &ps, const string &contig) const;
};
}
;
