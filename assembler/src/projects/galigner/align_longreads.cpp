//***************************************************************************
//* Copyright (c) 2015-2018 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>

#include "llvm/Support/YAMLParser.h"
#include "llvm/Support/YAMLTraits.h"
#include <cxxopts/cxxopts.hpp>

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
#include "pipeline/config_struct.hpp"

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/sequence_mapper.hpp"

#include "io/reads/io_helper.hpp"
#include "common/assembly_graph/core/coverage.hpp"

#include "modules/alignment/long_read_mapper.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "io/reads/multifile_reader.hpp"
#include "io/graph/gfa_reader.hpp"


#include "mapping_printer.hpp"

#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"


void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}


namespace graph_aligner {

struct GAlignerConfig {
    // general
    int K;
    string path_to_graphfile;
    string path_to_sequences;
    alignment::BWAIndex::AlignmentMode data_type; // pacbio, nanopore, 16S
    string output_format; // default: tsv and gpa without CIGAR

    //path construction
    debruijn_graph::config::pacbio_processor pb;
    GapClosingConfig gap_cfg;
};

}

namespace llvm {
namespace yaml {

template<> struct MappingTraits<debruijn_graph::config::pacbio_processor> {
    static void mapping(IO& io, debruijn_graph::config::pacbio_processor& cfg) {

        io.mapRequired("bwa_length_cutoff", cfg.bwa_length_cutoff);
        io.mapRequired("internal_length_cutoff", cfg.internal_length_cutoff);
        io.mapRequired("compression_cutoff", cfg.compression_cutoff);
        io.mapRequired("path_limit_stretching", cfg.path_limit_stretching);
        io.mapRequired("path_limit_pressing", cfg.path_limit_pressing);
        io.mapRequired("max_path_in_dijkstra", cfg.max_path_in_dijkstra);
        io.mapRequired("max_vertex_in_dijkstra", cfg.max_vertex_in_dijkstra);
        io.mapRequired("long_seq_limit", cfg.long_seq_limit);
        io.mapRequired("pacbio_min_gap_quantity", cfg.pacbio_min_gap_quantity);
        io.mapRequired("contigs_min_gap_quantity", cfg.contigs_min_gap_quantity);
        io.mapRequired("max_contigs_gap_length", cfg.max_contigs_gap_length);
    }
};

template<> struct MappingTraits<graph_aligner::GapClosingConfig> {
    static void mapping(IO& io, graph_aligner::GapClosingConfig& cfg) {
        io.mapRequired("max_vertex_in_gap", cfg.max_vertex_in_gap);
        io.mapRequired("queue_limit", cfg.queue_limit);
        io.mapRequired("iteration_limit", cfg.iteration_limit);
        io.mapRequired("find_shortest_path", cfg.find_shortest_path);
        io.mapRequired("restore_mapping", cfg.restore_mapping);
        io.mapRequired("penalty_interval", cfg.penalty_interval);
        io.mapRequired("max_ed_proportion", cfg.max_ed_proportion);
        io.mapRequired("ed_lower_bound", cfg.ed_lower_bound);
        io.mapRequired("ed_upper_bound", cfg.ed_upper_bound);
        io.mapRequired("max_restorable_end_length", cfg.max_restorable_end_length);
    }
};


template <>
struct ScalarEnumerationTraits<alignment::BWAIndex::AlignmentMode> {
    static void enumeration(IO &io, alignment::BWAIndex::AlignmentMode &value) {
        io.enumCase(value, "pacbio",  alignment::BWAIndex::AlignmentMode::PacBio);
        io.enumCase(value, "nanopore", alignment::BWAIndex::AlignmentMode::Ont2D);
        io.enumCase(value, "16S", alignment::BWAIndex::AlignmentMode::Rna16S);
    }
};

template<> struct MappingTraits<graph_aligner::GAlignerConfig> {
    static void mapping(IO& io, graph_aligner::GAlignerConfig& cfg) {
        io.mapRequired("k", cfg.K);
        io.mapRequired("path_to_graphfile", cfg.path_to_graphfile);
        io.mapRequired("path_to_sequences", cfg.path_to_sequences);
        io.mapRequired("data_type", cfg.data_type);
        io.mapRequired("output_format", cfg.output_format);
        io.mapRequired("run_dijkstra", cfg.gap_cfg.run_dijkstra);
        io.mapRequired("restore_ends", cfg.gap_cfg.restore_ends);

        io.mapRequired("pb", cfg.pb);
        io.mapRequired("gap_closing", cfg.gap_cfg);

    }
};

}
}

namespace graph_aligner {

class LongReadsAligner {
private:
    const debruijn_graph::ConjugateDeBruijnGraph &g_;
    const pacbio::PacBioMappingIndex<Graph> pac_index_;
    MappingPrinterHub mapping_printer_hub_;

    int aligned_reads_;
    int processed_reads_;

public:
    LongReadsAligner(const debruijn_graph::ConjugateDeBruijnGraph &g,
                     const alignment::BWAIndex::AlignmentMode mode,
                     const debruijn_graph::config::pacbio_processor &pb,
                     const GapClosingConfig gap_cfg,
                     const string output_file,
                     const string formats):
        g_(g), pac_index_(g_, pb, mode, gap_cfg), mapping_printer_hub_(g_, output_file, formats) {
        aligned_reads_ = 0;
        processed_reads_ = 0;
    }

    void AlignRead(const io::SingleRead &read) {
        DEBUG("Read " << read.name() << ". Current Read")
        utils::perf_counter pc;
        auto current_read_mapping = pac_index_.GetReadAlignment(read);
        const auto& aligned_mappings = current_read_mapping.main_storage;

        if (aligned_mappings.size() > 0) {
            mapping_printer_hub_.SaveMapping(current_read_mapping, read);
            DEBUG("Read " << read.name() << " is aligned");
            #pragma omp critical(align_read)
            {
                aligned_reads_ ++;
            }
        }  else {
            DEBUG("Read " << read.name() << " wasn't aligned and length=" << read.sequence().size());
        }
        #pragma omp critical(align_read)
        {
            processed_reads_ ++;
        }
        DEBUG("Read " << read.name() << " read_time=" << pc.time())
        return;
    }

    void RunAligner(std::vector<io::SingleRead> &wrappedreads, int threads) {
        size_t step = 10;
        processed_reads_ = 0;
        aligned_reads_ = 0;
        #pragma omp parallel num_threads(threads)
        #pragma omp for
        for (size_t i = 0 ; i < wrappedreads.size(); ++i) {
            AlignRead(wrappedreads[i]);
            #pragma omp critical(aligner)
            {
                if (processed_reads_ * 100 / wrappedreads.size() >= step) {
                    INFO("Processed \% reads: " << processed_reads_ * 100 / wrappedreads.size() <<
                         "\%, Aligned reads: " << aligned_reads_ * 100 / processed_reads_ <<
                         "\% (" << aligned_reads_ << " out of " << processed_reads_ << ")")
                    step += 10;
                }
            }
        }
    }
};

const debruijn_graph::ConjugateDeBruijnGraph& LoadGraph(const string &saves_path, const string &tmpdir, int K) {
    if (fs::extension(saves_path) == "gfa") {
        DEBUG("Load gfa")
        VERIFY_MSG(fs::is_regular_file(saves_path), "GFA-file " + saves_path + " doesn't exist");
        static debruijn_graph::ConjugateDeBruijnGraph g(K);
        gfa::GFAReader gfa(saves_path);
        DEBUG("Segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
        gfa.to_graph(g, true);
        return g;
    } else {
        DEBUG("Load from saves")
        static debruijn_graph::conj_graph_pack gp(K, tmpdir, 0);
        debruijn_graph::graphio::ScanGraphPack(saves_path, gp);
        return gp.g;
    }
}

void Launch(GAlignerConfig &cfg, const string output_file, int threads) {
    string tmpdir = fs::make_temp_dir(fs::current_dir(), "tmp");
    const debruijn_graph::ConjugateDeBruijnGraph &g = LoadGraph(cfg.path_to_graphfile, tmpdir, cfg.K);
    INFO("Loaded graph with " << g.size() << " vertices");

    LongReadsAligner aligner(g, cfg.data_type, cfg.pb, cfg.gap_cfg, output_file, cfg.output_format);
    INFO("LongReadsAligner created");

    io::ReadStreamList<io::SingleRead> streams;
    streams.push_back(make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(cfg.path_to_sequences)));
    io::SingleStreamPtr read_stream = io::MultifileWrap(streams);
    size_t n = 0;
    size_t buffer_no = 0;
    size_t read_buffer_size = 50000;
    while (!read_stream->eof()) {
        std::vector<io::SingleRead> read_buffer;
        read_buffer.reserve(read_buffer_size);
        io::SingleRead read;
        for (size_t buf_size = 0; buf_size < read_buffer_size && !read_stream->eof(); ++buf_size) {
            *read_stream >> read;
            read_buffer.push_back(std::move(read));
        }
        INFO("Prepared batch " << buffer_no << " of " << read_buffer.size() << " reads.");
        aligner.RunAligner(read_buffer, threads);
        ++buffer_no;
        n += read_buffer.size();
        INFO("Processed " << n << " reads");
    }
    INFO("Finished")
    fs::remove_dir(tmpdir);
}
} // namespace graph_aligner

int main(int argc, char **argv) {

    unsigned nthreads;
    std::string cfg, output_file;

    cxxopts::Options options(argv[0], " <YAML-config sequences and graph description> - Tool for sequence alignment on graph");
    options.add_options()
    ("o,outfile", "Output file prefix", cxxopts::value<std::string>(output_file)->default_value("./galigner_output"), "prefix")
    ("t,threads", "# of threads to use", cxxopts::value<unsigned>(nthreads)->default_value(std::to_string(min(omp_get_max_threads(), 16))), "num")
    ("h,help", "Print help");

    options.add_options("Input")
    ("positional", "", cxxopts::value<std::string>(cfg));

    options.parse_positional("positional");
    options.parse(argc, argv);
    if (options.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    if (!options.count("positional")) {
        std::cerr << "ERROR: No input YAML was specified" << std::endl << std::endl;
        std::cout << options.help() << std::endl;
        exit(-1);
    }

    create_console_logger();
    auto buf = llvm::MemoryBuffer::getFile(cfg);
    VERIFY_MSG(buf, "Failed to load config file " + cfg);
    llvm::yaml::Input yin(*buf.get());
    graph_aligner::GAlignerConfig config;
    yin >> config;
    omp_set_num_threads(nthreads);

    graph_aligner::Launch(config, output_file, nthreads);
    return 0;
}
