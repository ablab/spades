//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef COMMON_IO_BAMPARSER_HPP
#define COMMON_IO_BAMPARSER_HPP

#include "reads/single_read.hpp"
#include "io/reads/parser.hpp"
#include "sequence/quality.hpp"
#include "sequence/nucl.hpp"
#include "utils/verify.hpp"

#include "bamtools/api/BamReader.h"

#include <string>

namespace io {

class BAMParser: public Parser {
public:
    BAMParser(const std::string& filename, OffsetType offset_type = PhredOffset)
            : Parser(filename, offset_type) {
        open();
    }

    ~BAMParser() {
        close();
    }

    BAMParser& operator>>(SingleRead& read) {
        if (!is_open_ || eof_)
            return *this;

        read = SingleRead(seq_.Name, seq_.QueryBases, seq_.Qualities, offset_type_);
        eof_ = (false == reader_.GetNextAlignment(seq_));

        return *this;
    }

    void close() {
        reader_.Close();
        is_open_ = false;
        eof_ = true;
    }

private:
    BamTools::BamReader reader_;
    BamTools::BamAlignment seq_;

    void open() {
        reader_.Open(filename_);
        is_open_ = true;

        eof_ = (false == reader_.GetNextAlignment(seq_));
    }

    BAMParser(const BAMParser& parser);
    void operator=(const BAMParser& parser);
};

}

#endif /* COMMON_IO_FASTAFASTQGZPARSER_HPP */
