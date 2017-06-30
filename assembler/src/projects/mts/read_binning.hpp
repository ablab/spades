//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "annotation.hpp"
#include "io/reads/io_helper.hpp"
#include "gzstream/gzstream.h"

namespace io {

template<typename Stream>
class OSingleReadStream {
    Stream os_;

public:
    OSingleReadStream(const std::string& fn) {
        os_.open(fn.c_str());
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

template<typename Stream>
class OPairedReadStream {
    OSingleReadStream<Stream> l_os_;
    OSingleReadStream<Stream> r_os_;

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
    const EdgeAnnotation& edge_annotation_;
    std::string out_root_;
    std::string sample_name_;
    shared_ptr<SequenceMapper<Graph>> mapper_;
    std::set<std::string> bins_of_interest_;

    typedef io::OPairedReadStream<ogzstream> Stream;
    map<bin_id, std::shared_ptr<Stream>> out_streams_;

    set<bin_id> RelevantBins(const io::SingleRead& r) const;

    void Init(bin_id bin);

public:
    ContigBinner(const conj_graph_pack& gp,
                 const EdgeAnnotation& edge_annotation,
                 const std::string& out_root,
                 const std::string& sample_name,
                 const std::vector<std::string>& bins_of_interest = {}) :
                     gp_(gp),
                     edge_annotation_(edge_annotation),
                     out_root_(out_root),
                     sample_name_(sample_name),
                     mapper_(MapperInstance(gp)),
                     bins_of_interest_(bins_of_interest.begin(), bins_of_interest.end()) {
    }

    void Run(io::PairedStream& paired_reads);

    void close() {
        out_streams_.clear();
    }
};

int BinReads(const conj_graph_pack& gp, const std::string& out_root,
             const std::string& sample,
             const std::string& left_reads, const std::string& right_reads,
             const EdgeAnnotation& edge_annotation,
             const vector<string>& bins_of_interest);

}
