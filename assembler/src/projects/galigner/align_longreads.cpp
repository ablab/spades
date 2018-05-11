//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/standard_base.hpp"
#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
#include "pipeline/config_struct.hpp"

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/sequence_mapper.hpp"

#include "io/reads/io_helper.hpp"
#include "common/assembly_graph/core/coverage.hpp"

#include "modules/alignment/pacbio/pac_index.hpp"
#include "modules/alignment/long_read_mapper.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "io/reads/multifile_reader.hpp"
#include "mapping_printer.hpp"

#include "io/graph/gfa_reader.hpp"

#include "llvm/Support/YAMLParser.h"
#include "llvm/Support/YAMLTraits.h"
#include <cxxopts/cxxopts.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

struct GAlignerConfig {
    // general
    int K;
    string path_to_graphfile;
    string path_to_sequences; 
    string workdir;
    string data_type; // pacbio, nanopore, 16S
    
    string output_format; // default: tsv and gpa without CIGAR
    
    //path construction
    debruijn_graph::config::pacbio_processor pb;
    pacbio::GapClosingConfig gap_cfg;
};

}


namespace llvm { namespace yaml {

template<> struct MappingTraits<debruijn_graph::GAlignerConfig> {
    static void mapping(IO& io, debruijn_graph::GAlignerConfig& cfg) {
        io.mapRequired("k", cfg.K);
        io.mapRequired("path_to_graphfile", cfg.path_to_graphfile);
        io.mapRequired("path_to_sequences", cfg.path_to_sequences);
        io.mapRequired("workdir", cfg.workdir);
        io.mapRequired("data_type", cfg.data_type);
        io.mapRequired("output_format", cfg.output_format);
        io.mapRequired("run_dijkstra", cfg.gap_cfg.run_dijkstra);
        io.mapRequired("restore_ends", cfg.gap_cfg.restore_ends);

        io.mapRequired("pb.bwa_length_cutoff", cfg.pb.bwa_length_cutoff);
        io.mapRequired("pb.compression_cutoff", cfg.pb.compression_cutoff);
        io.mapRequired("pb.path_limit_stretching", cfg.pb.path_limit_stretching);
        io.mapRequired("pb.path_limit_pressing", cfg.pb.path_limit_pressing);
        io.mapRequired("pb.max_path_in_dijkstra", cfg.pb.max_path_in_dijkstra);
        io.mapRequired("pb.max_vertex_in_dijkstra", cfg.pb.max_vertex_in_dijkstra);
        io.mapRequired("pb.long_seq_limit", cfg.pb.long_seq_limit);
        io.mapRequired("pb.pacbio_min_gap_quantity", cfg.pb.pacbio_min_gap_quantity);
        io.mapRequired("pb.contigs_min_gap_quantity", cfg.pb.contigs_min_gap_quantity);
        io.mapRequired("pb.max_contigs_gap_length", cfg.pb.max_contigs_gap_length);

        io.mapRequired("gap_dijkstra.max_vertex_in_gap", cfg.gap_cfg.max_vertex_in_gap);
        io.mapRequired("gap_dijkstra.queue_limit", cfg.gap_cfg.queue_limit);
        io.mapRequired("gap_dijkstra.iteration_limit", cfg.gap_cfg.iteration_limit);
        io.mapRequired("gap_dijkstra.find_shortest_path", cfg.gap_cfg.find_shortest_path);
        io.mapRequired("gap_dijkstra.restore_mapping", cfg.gap_cfg.restore_mapping);
        io.mapRequired("gap_dijkstra.penalty_interval", cfg.gap_cfg.penalty_interval);
        
    }
};

} }

namespace debruijn_graph {

class BWASeedsAligner {
private:    
    const ConjugateDeBruijnGraph &g_;
    const pacbio::PacBioMappingIndex<Graph> pac_index_;
    MappingPrinterHub mapping_printer_hub_;

    int aligned_reads;
    int processed_reads;

public:
    BWASeedsAligner(const ConjugateDeBruijnGraph &g, 
                 const alignment::BWAIndex::AlignmentMode mode,
                 const debruijn_graph::config::pacbio_processor &pb,
                 const pacbio::GapClosingConfig gap_cfg,
                 const string output_file,
                 const string formats):
      g_(g),pac_index_(g_, pb, mode, gap_cfg), mapping_printer_hub_(g_, output_file, formats){
        aligned_reads = 0;
        processed_reads = 0;
      }

    void AlignRead(const io::SingleRead &read){
        DEBUG("Read " << read.name() <<". Current Read")
        auto current_read_mapping = pac_index_.GetReadAlignment(read);
        const auto& aligned_mappings = current_read_mapping.mapping_paths;
    
        if (aligned_mappings.size() > 0){
            mapping_printer_hub_.SaveMapping(aligned_mappings, read);
            DEBUG("Read " << read.name() <<" is aligned");
#pragma omp critical(align_read)
            {
                aligned_reads ++;
            }
        }  else {
            DEBUG("Read " << read.name() << " wasn't aligned and length=" << read.sequence().size());
        }
#pragma omp critical(align_read)
        {
            processed_reads ++;
        }
        return;
    } 

    void RunAligner(std::vector<io::SingleRead> &wrappedreads, int threads) {
        size_t step = 10;
#pragma omp parallel num_threads(threads)
#pragma omp for
        for (size_t i =0 ; i < wrappedreads.size(); ++i) {
            AlignRead(wrappedreads[i]);
#pragma omp critical(aligner)
            {
                if (processed_reads*100/wrappedreads.size() >= step) {
                    INFO("Processed \% reads: " << processed_reads*100/wrappedreads.size() << 
                         "\%, Aligned reads: " << aligned_reads*100/processed_reads << 
                         "\% (" << aligned_reads << " out of " << processed_reads << ")")
                    step += 10;
                }
            } 
        }
    } 
};

const ConjugateDeBruijnGraph& LoadGraph(const string saves_path, const string workdir, int K){
        if (saves_path.find(".gfa") != std::string::npos) {
            DEBUG("Load gfa")
            struct stat sb;
            VERIFY_MSG(stat(saves_path.c_str(), &sb) == 0 and S_ISREG(sb.st_mode), "GFA-file " + saves_path + " doesn't exist");   
            static ConjugateDeBruijnGraph g(K);
            gfa::GFAReader gfa(saves_path);
            DEBUG("Segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
            gfa.to_graph(g, true);
            return g;
        } else {
            DEBUG("Load from saves")
            struct stat sb;
            if (stat(workdir.c_str(), &sb) != 0 or !S_ISDIR(sb.st_mode))
            {
                int status = mkdir(workdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                VERIFY_MSG(status == 0, "Failed to create workdir " + workdir);   
            }
            static conj_graph_pack gp(K, workdir, 0);
            graphio::ScanGraphPack(saves_path, gp);
            return gp.g;
        }
}

void Launch(debruijn_graph::GAlignerConfig &cfg, const string output_file, int threads) {
    const ConjugateDeBruijnGraph &g = LoadGraph(cfg.path_to_graphfile, cfg.workdir, cfg.K);
    INFO("Loaded graph with " << g.size() << " vertices");
    io::ReadStreamList<io::SingleRead> streams;
    streams.push_back(make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(cfg.path_to_sequences)));
    
    // TODO - wrapped?
    io::SingleStreamPtr sstream = io::MultifileWrap(streams); 
    std::vector<io::SingleRead> wrappedreads; 
    while (!sstream->eof()) {
        io::SingleRead read;
        *sstream >> read;
        wrappedreads.push_back(std::move(read));
    }
    INFO("Loaded sequences from " << cfg.path_to_sequences);
    
    alignment::BWAIndex::AlignmentMode mode;
    if (cfg.data_type == "pacbio"){
        mode = alignment::BWAIndex::AlignmentMode::PacBio;
    } else if (cfg.data_type == "nanopore"){
        mode = alignment::BWAIndex::AlignmentMode::Ont2D;
    } else if (cfg.data_type == "16S"){
        mode = alignment::BWAIndex::AlignmentMode::Rna16S;
    } else {
        INFO("No appropriate mode")
        exit(-1);
    }
    
    BWASeedsAligner aligner(g, mode, cfg.pb, cfg.gap_cfg, output_file, cfg.output_format); 
    INFO("BWASeedsAligner created");

    INFO("Processing " << wrappedreads.size() << " reads..")
    aligner.RunAligner(wrappedreads, threads);
    INFO("Finished")
}
}

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
    debruijn_graph::GAlignerConfig config;
    yin >> config;
    omp_set_num_threads(nthreads);

    debruijn_graph::Launch(config, output_file, nthreads);
    return 0;
}
