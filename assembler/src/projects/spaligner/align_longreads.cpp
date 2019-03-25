//***************************************************************************
//* Copyright (c) 2018-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/utils/edge_namer.hpp"
#include "io/binary/graph.hpp"
#include "io/reads/io_helper.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "io/reads/multifile_reader.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/graph/gfa_writer.hpp"
#include "assembly_graph/core/graph.hpp"
#include "utils/logger/log_writers.hpp"
#include "modules/alignment/pacbio/g_aligner.hpp"

#include "mapping_printer.hpp"

#include "llvm/Support/YAMLParser.h"
#include "llvm/Support/YAMLTraits.h"

#include <iostream>
#include <fstream>
#include <cxxopts/cxxopts.hpp>

using namespace std;

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace llvm {
namespace yaml {

template<> struct MappingTraits<debruijn_graph::config::pacbio_processor> {
    static void mapping(IO& io, debruijn_graph::config::pacbio_processor& cfg) {
        io.mapRequired("internal_length_cutoff", cfg.internal_length_cutoff);
        io.mapRequired("compression_cutoff", cfg.compression_cutoff);
        io.mapRequired("path_limit_stretching", cfg.path_limit_stretching);
        io.mapRequired("path_limit_pressing", cfg.path_limit_pressing);
        io.mapRequired("max_path_in_dijkstra", cfg.max_path_in_dijkstra);
        io.mapRequired("max_vertex_in_dijkstra", cfg.max_vertex_in_dijkstra);
        io.mapRequired("long_seq_limit", cfg.long_seq_limit);
        io.mapRequired("pacbio_min_gap_quantity", cfg.pacbio_min_gap_quantity);
        io.mapRequired("contigs_min_gap_quantity", cfg.contigs_min_gap_quantity);
    }
};

template<> struct MappingTraits<sensitive_aligner::GapClosingConfig> {
    static void mapping(IO& io, sensitive_aligner::GapClosingConfig& cfg) {
        io.mapRequired("queue_limit", cfg.queue_limit);
        io.mapRequired("iteration_limit", cfg.iteration_limit);
        io.mapRequired("updates_limit", cfg.updates_limit);
        io.mapRequired("find_shortest_path", cfg.find_shortest_path);
        io.mapRequired("restore_mapping", cfg.restore_mapping);
        io.mapRequired("penalty_ratio", cfg.penalty_ratio);
        io.mapRequired("max_ed_proportion", cfg.max_ed_proportion);
        io.mapRequired("ed_lower_bound", cfg.ed_lower_bound);
        io.mapRequired("ed_upper_bound", cfg.ed_upper_bound);
        io.mapRequired("max_gs_states", cfg.max_gs_states);
    }
};

template<> struct MappingTraits<sensitive_aligner::EndsClosingConfig> {
    static void mapping(IO& io, sensitive_aligner::EndsClosingConfig& cfg) {
        io.mapRequired("queue_limit", cfg.queue_limit);
        io.mapRequired("iteration_limit", cfg.iteration_limit);
        io.mapRequired("updates_limit", cfg.updates_limit);
        io.mapRequired("find_shortest_path", cfg.find_shortest_path);
        io.mapRequired("restore_mapping", cfg.restore_mapping);
        io.mapRequired("penalty_ratio", cfg.penalty_ratio);
        io.mapRequired("max_ed_proportion", cfg.max_ed_proportion);
        io.mapRequired("ed_lower_bound", cfg.ed_lower_bound);
        io.mapRequired("ed_upper_bound", cfg.ed_upper_bound);
        io.mapRequired("max_restorable_length", cfg.max_restorable_length);
    }
};

template<> struct MappingTraits<sensitive_aligner::GAlignerConfig> {
    static void mapping(IO& io, sensitive_aligner::GAlignerConfig& cfg) {
        io.mapRequired("output_format", cfg.output_format);
        io.mapRequired("run_dijkstra", cfg.gap_cfg.run_dijkstra);
        io.mapRequired("restore_ends", cfg.restore_ends);

        io.mapRequired("pb", cfg.pb);
        io.mapRequired("gap_closing", cfg.gap_cfg);
        io.mapRequired("ends_recovering", cfg.ends_cfg);
    }
};

}
}

namespace sensitive_aligner {

class LongReadsAligner {
  public:
    LongReadsAligner(const debruijn_graph::ConjugateDeBruijnGraph &g,
                     const io::CanonicalEdgeHelper<debruijn_graph::Graph> &edge_namer,
                     const GAlignerConfig &cfg,
                     const string output_file,
                     const int threads)
        :g_(g),
         cfg_(cfg),
         galigner_(g_, cfg),
         threads_(threads),
         mapping_printer_hub_(g_, edge_namer, output_file, cfg.output_format) {
        aligned_reads_ = 0;
        processed_reads_ = 0;
    }

    void RunAligner() {
        io::ReadStreamList<io::SingleRead> streams;
        streams.push_back(make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(cfg_.path_to_sequences)));
        io::SingleStreamPtr read_stream = io::MultifileWrap(streams);
        size_t n = 0;
        size_t buffer_no = 0;
        while (!read_stream->eof()) {
            vector<io::SingleRead> read_buffer;
            read_buffer.reserve(read_buffer_size);
            io::SingleRead read;
            for (size_t buf_size = 0; buf_size < read_buffer_size && !read_stream->eof(); ++buf_size) {
                *read_stream >> read;
                read_buffer.push_back(move(read));
            }
            INFO("Prepared batch " << buffer_no << " of " << read_buffer.size() << " reads.");
            AlignBatch(read_buffer);
            ++buffer_no;
            n += read_buffer.size();
            INFO("Processed " << n << " reads");
        }
    }

  private:

    OneReadMapping AlignRead(const io::SingleRead &read) const {
        DEBUG("Read " << read.name() << ". Current Read")
        utils::perf_counter pc;
        auto current_read_mapping = galigner_.GetReadAlignment(read);
        const auto& aligned_mappings = current_read_mapping.edge_paths;
        if (aligned_mappings.size() > 0) {
            DEBUG("Read " << read.name() << " is aligned");
        }  else {
            DEBUG("Read " << read.name() << " wasn't aligned and length=" << read.sequence().size());
        }
        DEBUG("Read " << read.name() << " read_time=" << pc.time())
        return current_read_mapping;
    }

    void AlignBatch(const vector<io::SingleRead> &reads) {
        size_t step = 10;
        processed_reads_ = 0;
        aligned_reads_ = 0;
        #pragma omp parallel for schedule(guided, 50) num_threads(threads_)
        for (size_t i = 0 ; i < reads.size(); ++i) {
            OneReadMapping res = AlignRead(reads[i]);
            if (res.edge_paths.size() > 0) {
                mapping_printer_hub_.SaveMapping(res, reads[i]);
            }
            #pragma omp critical(aligner)
            {
                if (res.edge_paths.size() > 0) {
                    aligned_reads_ ++;
                }
                processed_reads_ ++;
                if (processed_reads_ * 100 / reads.size() >= step) {
                    INFO("Processed \% reads: " << processed_reads_ * 100 / reads.size() <<
                         "\%, Aligned reads: " << aligned_reads_ * 100 / processed_reads_ <<
                         "\% (" << aligned_reads_ << " out of " << processed_reads_ << ")")
                    step += 10;
                }
            }
        }
    }

    const size_t read_buffer_size = 50000;

    const debruijn_graph::ConjugateDeBruijnGraph &g_;
    const GAlignerConfig &cfg_;
    const sensitive_aligner::GAligner galigner_;
    const int threads_;
    MappingPrinterHub mapping_printer_hub_;

    int aligned_reads_;
    int processed_reads_;

};

void LoadGraph(const string &saves_path, debruijn_graph::ConjugateDeBruijnGraph &g, io::IdMapper<std::string> &id_mapper) {
    if (fs::extension(saves_path) == ".gfa") {
        DEBUG("Load gfa");
        VERIFY_MSG(fs::is_regular_file(saves_path), "GFA-file " + saves_path + " doesn't exist");
        gfa::GFAReader gfa(saves_path);
        DEBUG("Segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
        gfa.to_graph(g, &id_mapper);
        return;
    } else {
        DEBUG("Load from saves");
        io::binary::Load(saves_path, g);
        return;
    }
}

void Launch(GAlignerConfig &cfg, const string output_file, int threads) {
    string tmpdir = fs::make_temp_dir(fs::current_dir(), "tmp");
    debruijn_graph::ConjugateDeBruijnGraph g(cfg.K);
    io::IdMapper<std::string> id_mapper;
    LoadGraph(cfg.path_to_graphfile, g, id_mapper);
    io::CanonicalEdgeHelper<debruijn_graph::Graph> edge_namer(g, io::MapNamingF<debruijn_graph::Graph>(id_mapper));
    INFO("Loaded graph with " << g.size() << " vertices");

    LongReadsAligner aligner(g, edge_namer, cfg, output_file, threads);
    INFO("LongReadsAligner created");

    INFO("Process reads from " << cfg.path_to_sequences);
    aligner.RunAligner();
    INFO("Finished")
    fs::remove_dir(tmpdir);
}
} // namespace sensitive_aligner

int main(int argc, char **argv) {

    unsigned nthreads;
    string cfg, output_file, seq_type;
    sensitive_aligner::GAlignerConfig config;

    cxxopts::Options options(argv[0], " <YAML-config sequences and graph description> - Tool for sequence alignment on graph");
    options.add_options()
    ("K,k-mer", "graph k-mer size (odd value)", cxxopts::value<int>(config.K))
    ("d,datatype", "type of sequences: nanopore, pacbio", cxxopts::value<string>(seq_type))
    ("s,seq", "path to fasta/fastq file with sequences", cxxopts::value<string>(config.path_to_sequences))
    ("g,graph", "path to gfa-file or SPAdes saves folder", cxxopts::value<string>(config.path_to_graphfile))
    ("o,outfile", "Output file prefix", cxxopts::value<string>(output_file)->default_value("./galigner_output"), "prefix")
    ("t,threads", "# of threads to use", cxxopts::value<unsigned>(nthreads)->default_value(to_string(min(omp_get_max_threads(), 16))), "threads")
    ("h,help", "Print help");

    options.add_options("Input")
    ("positional", "", cxxopts::value<string>(cfg));

    options.parse_positional("positional");
    options.parse(argc, argv);
    if (options.count("help")) {
        cout << options.help() << endl;
        exit(0);
    }
    if (!options.count("positional")) {
        cerr << "ERROR: No input YAML was specified" << endl << endl;
        cout << options.help() << endl;
        exit(-1);
    }

    if (seq_type == "nanopore")
        config.data_type = alignment::BWAIndex::AlignmentMode::Ont2D;
    else if (seq_type == "pacbio")
        config.data_type = alignment::BWAIndex::AlignmentMode::PacBio;
    else {
        cerr << "ERROR: Unsupported data type - nanopore, pacbio" << endl;
        cout << options.help() << endl;
        exit(-1);
    }

    create_console_logger();
    INFO("Loading config from " << cfg)
    auto buf = llvm::MemoryBuffer::getFile(cfg);
    VERIFY_MSG(buf, "Failed to load config file " + cfg);
    llvm::yaml::Input yin(*buf.get());
    yin >> config;
    omp_set_num_threads(nthreads);

    sensitive_aligner::Launch(config, output_file, nthreads);
    return 0;
}
