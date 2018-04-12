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
#include "io/gfa/gfa_reader.hpp"

#include <iostream>
#include <fstream>

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

class BWASeedsAligner {
private:    
    const ConjugateDeBruijnGraph &g_;
    const pacbio::PacBioMappingIndex<Graph> pac_index_;
    MappingPrinterHub mapping_printer_hub_;

public:
    BWASeedsAligner(const ConjugateDeBruijnGraph &g, 
                 const alignment::BWAIndex::AlignmentMode mode,
                 const config::debruijn_config::pacbio_processor &pb,
                 const bool use_dijkstra,
                 const string output_file,
                 const string formats):
      g_(g),pac_index_(g_, pb, mode, use_dijkstra), mapping_printer_hub_(g_, output_file, formats){
      }

    // bool IsCanonical(EdgeId e) const {
    //     return e <= g.conjugate(e);
    // }

    // EdgeId Canonical(EdgeId e) const {
    //     return IsCanonical(e) ? e : g.conjugate(e);
    // }

    bool AlignRead(const io::SingleRead &read){
        Sequence seq(read.sequence());
        INFO("Read " << read.name() <<". Current Read")
        auto current_read_mapping = pac_index_.GetReadAlignment(seq);
        const auto& aligned_mappings = current_read_mapping.main_storage;
    
        if (aligned_mappings.size() > 0){
            mapping_printer_hub_.SaveMapping(aligned_mappings, read);
            INFO("Read " << read.name() <<" is aligned")
            return true;
        }  else {
            INFO("Read " << read.name() << " wasn't aligned and length=" << read.sequence().size());
            return false;
        }
    } 

    void RunAligner(std::vector<io::SingleRead> &wrappedreads, int threads) {
#pragma omp parallel num_threads(threads)
#pragma omp for
        for (size_t i =0 ; i < wrappedreads.size(); ++i) {
            AlignRead(wrappedreads[i]);
        }
    } 
};

config::debruijn_config::pacbio_processor InitializePacBioProcessor() {
    config::debruijn_config::pacbio_processor pb;  
    pb.bwa_length_cutoff = 0; //500
    pb.compression_cutoff = 0.6; // 0.6
    pb.path_limit_stretching = 1.3; //1.3
    pb.path_limit_pressing = 0.7;//0.7
    pb.max_path_in_dijkstra = 15000; //15000
    pb.max_vertex_in_dijkstra = 2000; //2000
//gap_closer
    pb.long_seq_limit = 400; //400
    pb.pacbio_min_gap_quantity = 2; //2
    pb.contigs_min_gap_quantity = 1; //1
    pb.max_contigs_gap_length = 10000; // 10000
    return pb;
}

debruijn_graph::ga_config InitializeGaConfig() {
    debruijn_graph::ga_config aligner_config;
    aligner_config.run_dijkstra = true;
    aligner_config.restore_edges = true;
    aligner_config.find_shortest_path = true;

    return aligner_config;
}

const ConjugateDeBruijnGraph& LoadGraphFromSaves(const string &saves_path, int K){
        if (saves_path.find(".gfa") != std::string::npos) {
            INFO("Load gfa")
            static ConjugateDeBruijnGraph g(K);
            gfa::GFAReader gfa(saves_path);
            INFO("Segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
            gfa.to_graph(g, true);
            return g;
        } else {
            INFO("Load from saves")
            static conj_graph_pack gp(K, "tmp3", 0);
            graphio::ScanGraphPack(saves_path, gp);
            return gp.g;
        }
}

void Launch(size_t K, const string &saves_path, 
                      const string &sequence_fasta, 
                      const string &mapper_type, 
                      const string &gapclose_type, 
                      const string &output_file, int threads) {
    const ConjugateDeBruijnGraph &g = LoadGraphFromSaves(saves_path, K);
    INFO("Loaded graph with " << g.size() << " vertices");
    io::ReadStreamList<io::SingleRead> streams;
    streams.push_back(make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(sequence_fasta)));
    
    io::SingleStreamPtr sstream = io::MultifileWrap(streams); 
    std::vector<io::SingleRead> wrappedreads; 
    while (!sstream->eof()) {
        io::SingleRead read;
        *sstream >> read;
        wrappedreads.push_back(std::move(read));
    }
    INFO("Loaded sequences from " << sequence_fasta);
    
    alignment::BWAIndex::AlignmentMode mode;
    if (mapper_type == "pacbio"){
        mode = alignment::BWAIndex::AlignmentMode::PacBio;
    } else if (mapper_type == "nanopore"){
        mode = alignment::BWAIndex::AlignmentMode::Ont2D;
    } else if (mapper_type == "16S"){
        mode = alignment::BWAIndex::AlignmentMode::Rna16S;
    } else {
        INFO("No appropriate mode")
        mode = alignment::BWAIndex::AlignmentMode::PacBio;
        exit(-1);
    }

    debruijn_graph::ga_config aligner_config = InitializeGaConfig();

    config::debruijn_config::pacbio_processor pb = InitializePacBioProcessor();
    bool use_dijkstra = true;
    if (gapclose_type == "bf") {
        use_dijkstra = false;
    }
    BWASeedsAligner aligner(g, mode, pb, use_dijkstra, output_file, "tsv"); 
    INFO("BWASeedsAligner created");

    aligner.RunAligner(wrappedreads, threads);
}
}

int main(int argc, char **argv) {
    if (argc < 7) {
        cout << "Usage: longreads_aligner <K>"
             << " <saves path> <sequences file (fasta/fastq)> <mapper type: {pacbio, nanopore, 16S} > <gap close type: {bf, dijkstra(default)} > <ouput-prefix> {threads-num(16)}" << endl;
        exit(1);
    }
    int threads = 16;
    if (argc == 8) {
        threads = std::stoi(argv[7]);
    }
    omp_set_num_threads(threads);
    create_console_logger();
    size_t K = std::stoll(argv[1]);
    string saves_path = argv[2];
    INFO("Load graph from " << saves_path);
    string sequence_file = argv[3];
    INFO("Load sequences from " << sequence_file);
    string mapper_type = argv[4];
    INFO("Mapper type " << mapper_type);
    string gapclose_type = argv[5];
    INFO("Gap close type " << gapclose_type);
    string output_file = argv[6];
    debruijn_graph::Launch(K, saves_path, sequence_file, mapper_type, gapclose_type, output_file, threads);
    return 0;
}
