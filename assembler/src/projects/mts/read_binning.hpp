//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "annotation.hpp"
#include "io/reads_io/io_helper.hpp"

namespace io {

class OSingleReadStream {
    std::ofstream os_;

public:
    OSingleReadStream(const std::string& fn) :
        os_(fn) {
    }

    OSingleReadStream& operator<<(const SingleRead& read) {
        os_ << "@" << read.name() << std::endl;
        os_ << read.GetSequenceString() << std::endl;
        os_ << "+" << std::endl;
        os_ << read.GetPhredQualityString() << std::endl;
        return *this;
    }

    void close() {
        os_.close();
    }
};

class OPairedReadStream {
    OSingleReadStream l_os_;
    OSingleReadStream r_os_;

public:
    OPairedReadStream(const std::string& l_fn, const std::string& r_fn) :
        l_os_(l_fn), r_os_(r_fn) {
    }

    OPairedReadStream& operator<<(const PairedRead& read) {
        l_os_ << read.first();
        r_os_ << read.second();
        return *this;
    }

    void close() {
        l_os_.close();
        r_os_.close();
    }
};

}

namespace debruijn_graph {

class ContigBinner {
    const conj_graph_pack& gp_;
    EdgeAnnotation edge_annotation_;

    map<bin_id, std::shared_ptr<io::OPairedReadStream>> out_streams_;

public:
    ContigBinner(const conj_graph_pack& gp, const vector<bin_id>& bins_of_interest) :
                     gp_(gp),
                     edge_annotation_(gp, bins_of_interest) {
    }

    void Init(const string& output_root, const string& sample_name, io::SingleStream& contigs, AnnotationStream& annotation_stream);

    void Run(io::PairedStream& paired_reads);

    void close() {
        out_streams_.clear();
    }
};

}
