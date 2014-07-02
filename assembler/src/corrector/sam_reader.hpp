//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
//todo rename to reader
#pragma once


#include <samtools/sam.h>
#include "samtools/bam.h"
#include "read.hpp"
#include "io/ireader.hpp"

//namespace io {
/*class SamRead : public SamTools::BamAlignment {
  public:
    BamRead() {}

    BamRead(const BamTools::BamAlignment &other)
            : BamTools::BamAlignment(other) {}

    const std::string& name() const {
        return Name;
    }

    size_t size() const {
        return Length;
    }

    size_t nucl_count() const {
        return size();
    }

    const std::string& GetSequenceString() const {
        return QueryBases;
    }

    std::string GetPhredQualityString() const {
        return Qualities;
    }

    operator io::SingleRead() {
        // not including quality is intentional:
        // during read correction bases might be inserted/deleted,
        // and base qualities for them are not calculated
        return io::SingleRead(name(), GetSequenceString());
    }

    char operator[](size_t i) const {
        VERIFY(is_nucl(QueryBases[i]));
        return dignucl(QueryBases[i]);
    }
};
*/
class MappedSamStream: public io::ReadStream<SingleSamRead> {
  public:
    MappedSamStream(const std::string &filename)
            : filename_(filename) {
    	open();
    }

    virtual ~MappedSamStream() {}

    bool is_open() { return is_open_; }
    bool eof() { return eof_; }
    MappedSamStream& operator>>(SingleSamRead& read) {
        if (!is_open_ || eof_)
            return *this;

     //   read = seq_;
//        cerr <<"name: " <<  bam1_qname(seq_);
//        //cerr << " " << cigar(bam1_cigar(seq_));
//        cerr << " cigar: " << GetCigar(seq_);
//        cerr << " seq: " << GetSeq(seq_);
//        cerr << " qual: " << GetQual(seq_);
//        cerr << '\n';
        unsigned char* ls = seq_->data;
        read.data_ = seq_;
        int tmp = samread(reader_, seq_);
        eof_ = (0 >= tmp);
        return *this;
    }
    void ReadHeader(bam_header_t *bam_header){
//    	INFO(reader_->header->n_targets << " contigs in header ");
//    	INFO(reader_->header->target_name[0]);
    	bam_header = reader_->header;
    }

    string get_contig_name(int i){
    	VERIFY(i < reader_->header->n_targets);
    	return (reader_->header->target_name[i]);
    }
    void close() {
    	samclose(reader_);
    	is_open_ = false;
        eof_ = true;
    }

    void reset() {
        close();
        open();
    }

    io::ReadStreamStat get_stat() const { return io::ReadStreamStat(); }

  private:
    samfile_t *reader_;
    bam1_t *seq_ = new bam1_t;
    std::string filename_;
    bool is_open_;
    bool eof_;


    void open() {
        if ((reader_ = samopen(filename_.c_str(), "r", NULL)) == NULL)
		{
		   cerr << "Fail to open SAM/BAM file " << filename_ << endl;
		 //  exit(-1);
		}
        is_open_ = true;
        int tmp = samread(reader_, seq_);
        eof_ = (0 >= tmp);

    }

};
//}
