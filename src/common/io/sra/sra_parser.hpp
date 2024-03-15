//***************************************************************************
//* Copyright (c) 2023 SPAdes authors
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/file_read_flags.hpp"
#include "io/reads/single_read.hpp"
#include "io/reads/parser.hpp"

namespace ngs {
class ReadCollection;
class ReadIterator;
};

#include <string>

namespace io {

class SRAParser: public Parser {
public:
    SRAParser(const std::string& filename,
              FileReadFlags flags = FileReadFlags())
            : Parser(filename, flags) {
        open();
    }

    ~SRAParser() {
        close();
    }

    SRAParser& operator>>(SingleRead& read);
    void close();

    SRAParser(const SRAParser& parser) = delete;
    void operator=(const SRAParser& parser) = delete;

private:
    void open();
    void next();

    ngs::ReadCollection* run_;
    ngs::ReadIterator* it_;
};

}
