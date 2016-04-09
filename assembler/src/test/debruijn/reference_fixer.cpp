//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dev_support/standard_base.hpp"
#include "dev_support/simple_tools.hpp"
#include "dev_support/logger/log_writers.hpp"

#include "dev_support/path_helper.hpp"
#include "io/reads_io/file_reader.hpp"
#include "io/reads_io/wrapper_collection.hpp"
#include "io/reads_io/osequencestream.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cout << "Usage: reference_fixer <reference path> <output path>" << endl;
        return 1;
    }
    create_console_logger();
    string fn = argv[1];
    path::CheckFileExistenceFATAL(fn);
    string out_fn = argv[2];
    auto reader = make_shared<io::NonNuclCollapsingWrapper>(make_shared<io::FileReadStream>(fn));
    io::SingleRead read;
    std::stringstream ss;
    while (!reader->eof()) {
        (*reader) >> read;
        ss << read.GetSequenceString();
    }
    io::SingleRead concat("concat", ss.str());
    io::osequencestream out(out_fn);
    out << concat;
    return 0;
}
