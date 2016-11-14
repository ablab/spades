//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <io/sam/read.hpp>
#include <io/sam/sam_reader.hpp>

namespace sam_reader {

bool MappedSamStream::eof() const {
        return eof_;
}

bool MappedSamStream::is_open() const {
    return is_open_;
}

MappedSamStream& MappedSamStream::operator>>(SingleSamRead& read) {
    if (!is_open_ || eof_)
        return *this;
    read.set_data(seq_);
    int tmp = samread(reader_, seq_);
    eof_ = (0 >= tmp);
    return *this;
}

MappedSamStream& MappedSamStream::operator >>(PairedSamRead& read) {
    TRACE("starting process paired read");
    SingleSamRead r1;
    MappedSamStream::operator >>(r1);
    SingleSamRead r2;
    MappedSamStream::operator >>(r2);

    read = PairedSamRead(r1, r2);
    TRACE(r1.seq());
    TRACE(r2.seq());
    TRACE(r1.name());
    return *this;
}

const char* MappedSamStream::get_contig_name(int i) const {
    VERIFY(i < reader_->header->n_targets);
    return (reader_->header->target_name[i]);
}

void MappedSamStream::close() {
    samclose(reader_);
    is_open_ = false;
    eof_ = true;
    bam_destroy1(seq_);
}

void MappedSamStream::reset() {
    close();
    open();
}

void MappedSamStream::open() {
    if ((reader_ = samopen(filename_.c_str(), "r", NULL)) == NULL) {
        WARN("Fail to open SAM file " << filename_);
        is_open_ = false;
        eof_ = true;
    } else {
        is_open_ = true;
        int tmp = samread(reader_, seq_);
        eof_ = (0 >= tmp);
    }
}

}
