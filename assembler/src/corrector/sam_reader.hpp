//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
#pragma once


// WTF: Make sure your includes are in proper order
#include <samtools/sam.h>
#include "samtools/bam.h"
#include "read.hpp"
#include "io/ireader.hpp"

// WTF: EVERYWHERE: USE SPACES, NOT TABS! FIX ALL THE CODING STYLE PROBLEMS EVERYWHERE

namespace corrector {

class MappedSamStream: public io::ReadStream<SingleSamRead> {
  public:
    MappedSamStream(const std::string &filename)
            : filename_(filename) {
    	open();
    }

    virtual ~MappedSamStream() {}

    // WTF: Why these are not const?
    bool is_open();
    bool eof();
    MappedSamStream& operator>>(SingleSamRead& read);
    MappedSamStream& operator >> (PairedSamRead& read);
    bam_header_t* ReadHeader();
    string get_contig_name(int i);
    void close();
    void reset();
    io::ReadStreamStat get_stat() const;

  private:
    samfile_t *reader_;
    bam1_t *seq_ =bam_init1();
    std::string filename_;
    bool is_open_;
    bool eof_;


    void open();
};

};
