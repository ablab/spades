//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "bam_parser.hpp"

#include <bamtools/api/BamReader.h>

namespace io {

BAMParser& BAMParser::operator>>(SingleRead& read) {
    if (!is_open_ || eof_)
        return *this;

    if (flags_.use_name && flags_.use_quality)
        read = SingleRead(seq_.Name, "", seq_.QueryBases, seq_.Qualities, flags_.offset,
                          0, 0, flags_.validate);
    else if (flags_.use_name)
        read = SingleRead(seq_.Name, "", seq_.QueryBases,
                          0, 0, flags_.validate);
    else
        read = SingleRead(seq_.QueryBases,
                          0, 0, flags_.validate);

    eof_ = (false == reader_.GetNextAlignment(seq_));

    return *this;
}

void BAMParser::close() {
    reader_.Close();
    is_open_ = false;
    eof_ = true;
}

void BAMParser::open() {
    reader_.Open(filename_);
    is_open_ = true;

    eof_ = (false == reader_.GetNextAlignment(seq_));
}

}
