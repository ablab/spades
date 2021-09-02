//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sra_parser.hpp"

#include <ncbi-vdb/NGS.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>

namespace io {

void SRAParser::open() {
    run_ = new ngs::ReadCollection(ncbi::NGS::openReadCollection(filename_));
    it_ = new ngs::ReadIterator(run_->getReads(ngs::Read::all));
    is_open_ = true;
    eof_ = (false == it_->nextRead());
    while (!eof_ && !it_->nextFragment()) {
        eof_ = (false == it_->nextRead());
    }
}

SRAParser& SRAParser::operator>>(SingleRead& read) {
    if (!is_open_ || eof_)
        return *this;

    if (flags_.use_name && flags_.use_quality)
        read = SingleRead(it_->getReadName().toString(), "",
                          it_->getFragmentBases().toString(), it_->getFragmentQualities().toString(), flags_.offset,
                          0, 0, flags_.validate);
    else if (flags_.use_name)
        read = SingleRead(it_->getReadName().toString(), "",
                          it_->getFragmentBases().toString(),
                          0, 0, flags_.validate);
    else
        read = SingleRead(it_->getFragmentBases().toString(),
                          0, 0, flags_.validate);

    while (!eof_ && !it_->nextFragment()) {
        eof_ = (false == it_->nextRead());
    }
    
    return *this;
}

void SRAParser::close() {
    is_open_ = false;
    eof_ = true;
    delete it_;
    delete run_;
}

}
