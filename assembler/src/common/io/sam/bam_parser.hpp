//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/file_read_flags.hpp"
#include "io/reads/single_read.hpp"
#include "io/reads/parser.hpp"

#include "bamtools/api/BamReader.h"

#include <string>

namespace io {

class BAMParser: public Parser {
public:
    BAMParser(const std::string& filename,
              FileReadFlags flags = FileReadFlags())
            : Parser(filename, flags) {
        open();
    }

    ~BAMParser() {
        close();
    }

    BAMParser& operator>>(SingleRead& read);
    void close();

    BAMParser(const BAMParser& parser) = delete;
    void operator=(const BAMParser& parser) = delete;

private:
    BamTools::BamReader reader_;
    BamTools::BamAlignment seq_;

    void open();
};

}
