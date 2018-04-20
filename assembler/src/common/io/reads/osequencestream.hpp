//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "single_read.hpp"
#include "paired_read.hpp"
#include "header_naming.hpp"

#include <fstream>
#include <string>
#include <vector>

namespace io {

inline void WriteWrapped(const std::string &s, std::ostream &os, size_t max_width = 60) {
    size_t cur = 0;
    while (cur < s.size()) {
        os << s.substr(cur, max_width) << "\n";
        cur += max_width;
    }
}

class osequencestream {
protected:
    std::ofstream ofstream_;
    size_t id_;

    void write_str(const std::string& s) {
        WriteWrapped(s, ofstream_);
    }

    virtual void write_header(const std::string& s) {
        // Velvet format: NODE_1_length_24705_cov_358.255249
        ofstream_ << ">" << MakeContigId(id_++, s.size()) << std::endl;
    }

public:
    osequencestream(const std::string& filename):
            ofstream_(filename), id_(1) {
    }

    virtual ~osequencestream() {}

    osequencestream& operator<<(const std::string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

    osequencestream& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

};

class osequencestream_cov: public osequencestream {
protected:
    double coverage_;

    virtual void write_header(const std::string& s) {
        // Velvet format: NODE_1_length_24705_cov_358.255249
        ofstream_ << ">" << MakeContigId(id_++, s.size(), coverage_) << std::endl;
    }


public:
    osequencestream_cov(const std::string& filename)
        : osequencestream(filename), coverage_(0.) { }

    osequencestream_cov& operator<<(double coverage) {
        coverage_ = coverage;
        return *this;
    }

    using osequencestream::operator<<;

};

class OFastaReadStream {
    std::ofstream ofstream_;
public:
    typedef SingleRead ReadT;

    OFastaReadStream(const std::string& filename):
            ofstream_(filename) {
    }

    OFastaReadStream& operator<<(const SingleRead& read) {
        ofstream_ << ">" << read.name() << "\n";
        WriteWrapped(read.GetSequenceString(), ofstream_);
        return *this;
    }
};

class OFastqReadStream {
    std::ofstream os_;

public:
    typedef SingleRead ReadT;

    OFastqReadStream(const std::string& fn) :
            os_(fn) {
    }

    OFastqReadStream& operator<<(const SingleRead& read) {
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

template<class SingleReadStream>
class OPairedReadStream {
    SingleReadStream l_os_;
    SingleReadStream r_os_;

public:
    typedef PairedRead ReadT;

    OPairedReadStream(const std::string& l_fn,
                            const std::string& r_fn) :
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

using OFastaPairedStream = OPairedReadStream<OFastaReadStream>;
using OFastqPairedStream = OPairedReadStream<OFastqReadStream>;

}
