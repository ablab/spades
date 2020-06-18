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

    void write_header(const std::string& s) override {
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

class osequencestream_bgc: public osequencestream {
protected:
    unsigned cluster_;
    unsigned candidate_;
    unsigned domains_;

    virtual void write_header(const std::string& s) {
        // Velvet format: NODE_1_length_24705_cov_358.255249
        ofstream_ << ">" << AddClusterId(MakeContigId(id_++, s.size()), cluster_, candidate_, domains_) << std::endl;
    }


public:
    osequencestream_bgc(const std::string& filename)
            : osequencestream(filename), cluster_(0), candidate_(0) { }

    void SetCluster(unsigned cluster, unsigned candidate, unsigned domains) {
        cluster_ = cluster;
        candidate_ = candidate;
        domains_ = domains;
    }

    using osequencestream::operator<<;
};

struct FastaWriter {
    static void Write(std::ostream &stream, const SingleRead &read) {
        stream << ">" << read.name() << "\n";
        WriteWrapped(read.GetSequenceString(), stream);
    }
};

struct FastqWriter {
    static void Write(std::ostream &stream, const SingleRead &read) {
        stream << "@" << read.name() << std::endl
               << read.GetSequenceString() << std::endl
               << "+" << std::endl
               << read.GetPhredQualityString() << std::endl;
    }
};

template<typename Stream, typename Writer>
class OReadStream {
public:
    typedef SingleRead ReadT;

    OReadStream(const std::string &filename)
            : stream_(filename) {
    }

    OReadStream &operator<<(const SingleRead &read) {
        Writer::Write(stream_, read);
        return *this;
    }

    void close() {
        stream_.close();
    }

private:
    Stream stream_;
};

typedef OReadStream<std::ofstream, FastaWriter> OFastaReadStream;
typedef OReadStream<std::ofstream, FastqWriter> OFastqReadStream;

template<typename Stream, typename Writer>
class OPairedReadStream {
public:
    typedef PairedRead ReadT;

    OPairedReadStream(const std::string &left_filename, const std::string &right_filename)
            : left_stream_(left_filename), right_stream_(right_filename) {
    }

    OPairedReadStream &operator<<(const PairedRead &read) {
        Writer::Write(left_stream_, read.first());
        Writer::Write(right_stream_, read.second());
        return *this;
    }

    void close() {
        left_stream_.close();
        right_stream_.close();
    }

private:
    Stream left_stream_, right_stream_;
};

typedef OPairedReadStream<std::ofstream, FastaWriter> OFastaPairedStream;
typedef OPairedReadStream<std::ofstream, FastqWriter> OFastqPairedStream;

}
