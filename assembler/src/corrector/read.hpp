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
        //Re: no new->nothing to delete. data_ is freed in read destructor, *seq_ is freed in sam_reader, new_seq - right after this function
        bam1_t *new_seq = bam_dup1(seq_);
        data_ = *new_seq;
    }
    // WTF: This does not belong here
    // Re: Now it knows nothing about the contig, and depends on read and length parameter only. I'm still sure this should be here- bam_1 stuff should be localized.
    int CountPositions(unordered_map<size_t, position_description> &ps, const size_t &contig_length) const;
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
    // Re: same as for SingleSamRead
    int CountPositions(unordered_map<size_t, position_description> &ps, const size_t &contig_length) const;
};
}
;
