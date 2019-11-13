//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "version.hpp"

#include "io/reads/read_processor.hpp"
#include "io/reads/io_helper.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/ph_map/perfect_hash_map.hpp"
#include "utils/kmer_mph/kmer_index_builder.hpp"
#include "pipeline/library.hpp"

#include <clipp/clipp.h>
#include <string>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cerrno>

#include <sys/types.h>
#include <sys/stat.h>
#include <common/io/dataset_support/read_converter.hpp>
#include <common/io/reads/osequencestream.hpp>

using namespace std;
void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

namespace convert_bin_to_fasta {
struct Args {
    std::string prefix_file = "";
    std::string info_file = "";
    std::string output_file = "contigs.fasta";
};
}

void process_cmdline(int argc, char **argv, convert_bin_to_fasta::Args &args) {
    using namespace clipp;
    bool print_help = false;

    auto cli = (
        (option("--prefix") & value("path", args.prefix_file)) % "Prefix of .off and .seq file for contigs in binary format",
        (option("--info_file") & value("path", args.info_file)) % "Path to info file for contigs in binary format",
        (option("-o", "--output_file") & value("path", args.output_file)) % "Output file name",
        (option("-h", "--help").set(print_help)) % "Show help"
    );

    auto help_message = make_man_page(cli, argv[0])
        .prepend_section("DESCRIPTION",
                         "Convert contigs in binary format to fasta format.");

    auto result = parse(argc, argv, cli);
    if (!result || print_help) {
        std::cout << help_message;
        if (print_help) {
            exit(0);
        } else {
            exit(1);
        }
    }

    if (args.prefix_file == "" || args.info_file == "") {
        std::cerr << "ERROR: No binary file were specified (you should specify info file and prefix to .off and .seq)"
                  << std::endl << std::endl;
        std::cout << help_message << std::endl;
        exit(-1);
    }
}

io::ReadStreamList<io::SingleReadSeq> get_bin_stream(const convert_bin_to_fasta::Args& args) {
    io::SequencingLibraryT seq_lib;
    seq_lib.set_type(io::LibraryType::TrustedContigs);
    seq_lib.set_orientation(io::LibraryOrientation::Undefined);
    seq_lib.data().lib_index = size_t(-1);
    auto& bin_info = seq_lib.data().binary_reads_info;
    bin_info.single_read_prefix = args.prefix_file;
    bin_info.bin_reads_info_file = args.info_file;
    bin_info.binary_converted = true;
    bin_info.chunk_num = 1;

    io::ReadStreamList<io::SingleReadSeq> lib_streams = io::single_binary_readers(seq_lib, true, false);
    return lib_streams;
}

void read_from_stream(io::BinarySingleStream &stream, io::osequencestream& oss) {
    io::SingleReadSeq r;
    while (!stream.eof()) {
        stream >> r;
        oss << r.sequence();
    }
}

int main(int argc, char* argv[]) {
    try {
        convert_bin_to_fasta::Args args;
        process_cmdline(argc, argv, args);

        create_console_logger();

        START_BANNER("SPAdes converter from binary to fasta");

        INFO("Binary files: " << args.prefix_file << " " << args.info_file);
        INFO("Output file: " << args.output_file);

        auto streams = get_bin_stream(args);
        io::osequencestream oss(args.output_file);

        streams.reset();
        while (!streams.eof()) {
            for (unsigned i = 0; i < streams.size(); ++i) {
                read_from_stream(streams[i], oss);
            }
        }
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    }

    return 0;
}
