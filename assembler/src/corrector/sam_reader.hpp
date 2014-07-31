//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
#pragma once


#include <samtools/sam.h>
#include "samtools/bam.h"
#include "read.hpp"
#include "io/ireader.hpp"

namespace corrector {

class MappedSamStream: public io::ReadStream<SingleSamRead> {
  public:
    MappedSamStream(const std::string &filename)
            : filename_(filename) {
    	open();
    }

    virtual ~MappedSamStream() {}

    bool is_open();
    bool eof();

    MappedSamStream& operator>>(SingleSamRead& read);/* {
        if (!is_open_ || eof_)
            return *this;
        read.set_data(seq_);
        int tmp = samread(reader_, seq_);
        eof_ = (0 >= tmp);
        return *this;
    }*/

    MappedSamStream& operator >> (PairedSamRead& read);/*{
    	TRACE("starting process paired read");
    	SingleSamRead r1;
    	MappedSamStream::operator >> (r1);
    	SingleSamRead r2;
    	MappedSamStream::operator >> (r2);
    	TRACE(r1.GetSeq());
    	TRACE(r2.GetSeq());
    	TRACE(r1.GetName());
    	VERIFY_MSG (r1.GetName() == r2.GetName(), r1.GetName() + " " + r2.GetName());
    	read.pair(r1,r2);
        return *this;
    }
*/
    bam_header_t* ReadHeader();/*{
    	return reader_->header;
    }*/

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
//}
};
