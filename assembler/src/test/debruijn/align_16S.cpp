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
#include "modules/alignment/short_read_mapper.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "io/reads/multifile_reader.hpp"

#include <iostream>
#include <fstream>

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

struct ReadMapping {

    void ExtractPositions() {
        start_pos_ = -1;
        end_pos_ = -1;
        string name = read_.name();
        int open_bracket_ind = -1;
        int close_bracket_ind = -1;
        int comma_ind = -1;
        for (int i = 0; i < name.size(); ++ i) {
            if (name[i] == '[') {
                open_bracket_ind = i;
            }
            if (name[i] == ']') {
                close_bracket_ind = i;
            }
            if (name[i] == ',') {
                comma_ind = i;
            }
        }
        start_pos_ = std::stoi(name.substr(open_bracket_ind + 1, comma_ind - open_bracket_ind));
        end_pos_ = std::stoi(name.substr(comma_ind + 1, close_bracket_ind - comma_ind));
        INFO("Name=" << name << " start_pos_=" << start_pos_ << " end_pos_=" << end_pos_);
    }

    ReadMapping(const io::SingleRead &read, std::vector<int> &v, std::vector<EdgeId> &e, std::vector<MappingRange> &range) 
        : read_(read), v_(v), e_(e), range_(range) {
        ExtractPositions();
    }

    io::SingleRead read_;
    std::vector<int> v_;
    std::vector<EdgeId> e_;
    std::vector<MappingRange> range_;
    int start_pos_;
    int end_pos_;
};


class PacBioAligner {
private:    
    const conj_graph_pack &gp_;
    const pacbio::PacBioMappingIndex<Graph> pac_index_;
    const string &output_file_;

public:
    PacBioAligner(const conj_graph_pack &gp, 
                 const alignment::BWAIndex::AlignmentMode mode,
                 const config::debruijn_config::pacbio_processor &pb,
                 const string &output_file):
      gp_(gp),pac_index_(gp.g, "./tmp4", pb, mode), output_file_(output_file){}

    
    bool IsCanonical(EdgeId e) const {
        return e <= gp_.g.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : gp_.g.conjugate(e);
    }

    int EditDistance(const string &a, const string &b, int &start_pos, int &end_pos) {
        int a_len = (int) a.length();
        int b_len = (int) b.length();
        VERIFY(a_len > 0);
        VERIFY(b_len > 0);
        edlib::EdlibEqualityPair additionalEqualities[36] = {{'U', 'T'}
                                                , {'R', 'A'}, {'R', 'G'}
                                                , {'Y', 'C'}, {'Y', 'T'}, {'Y', 'U'}
                                                , {'K', 'G'}, {'K', 'T'}, {'K', 'U'}
                                                , {'M', 'A'}, {'M', 'C'}
                                                , {'S', 'C'}, {'S', 'G'}
                                                , {'W', 'A'}, {'W', 'T'}, {'W', 'U'}
                                                , {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'B', 'U'}
                                                , {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'D', 'U'}
                                                , {'H', 'A'}, {'H', 'C'}, {'H', 'T'}, {'H', 'U'}
                                                , {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
                                                , {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'N', 'U'} };
        edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len, b.c_str(), b_len
                                                       , edlib::edlibNewAlignConfig(a_len, edlib::EDLIB_MODE_HW, edlib::EDLIB_TASK_LOC,
                                                                             additionalEqualities, 36));
        int score = -1;
        if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
            score = result.editDistance;
            start_pos = result.startLocations[0];
            end_pos = result.endLocations[0];
        }
        edlib::edlibFreeAlignResult(result);
        return score;
    }

    void AlignPrimer(const io::SingleRead &read, std::vector<int> &v, std::vector<EdgeId> &e, std::vector<MappingRange> &range) {
        for (int i = 0; i < v.size(); ++ i) {
            v[i] = -1;
            e[i] = EdgeId();
            range.push_back(MappingRange(-1, -1, -1, -1));
        }
        for (auto it = gp_.g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId eid = *it;
            std::string edge_str = gp_.g.EdgeNucls(eid).str();
            int start_pos = -1;
            int end_pos = -1;
            int dist = EditDistance(read.sequence().str(), edge_str, start_pos, end_pos);
            if (dist == -1) continue;
            int mx_ind = 0;
            for (int i = 0; i < v.size(); ++ i) {
                if (v[i] > v[mx_ind] || v[i] == -1) {
                    mx_ind = i;
                    if (v[i] == -1) break;
                }
            }
            if (v[mx_ind] == -1 || v[mx_ind] > dist) {
                v[mx_ind] = dist;
                e[mx_ind] = eid;
                range[mx_ind] = MappingRange(Range(0, read.size()), Range(start_pos, end_pos + 1));
            }
        }
        std::string ans = "";
        int mn_ind = 0;
        int mx_ind = 0;
        for (int i = 0; i < v.size(); ++ i) {
            //ans += " edge_id=" + std::to_string(e[i].int_id()) + " dist=" +  std::to_string(v[i]) + ";";
            if (v[i] > v[mx_ind]) {
                mx_ind = i;
            }
            if (v[i] < v[mn_ind]) {
                mn_ind = i;
            }
        }
        ans = " max: edge_id=" + std::to_string(e[mx_ind].int_id()) + " dist=" +  std::to_string(v[mx_ind]) + ";"
                + " min: edge_id=" + std::to_string(e[mn_ind].int_id()) + " dist=" +  std::to_string(v[mn_ind]);
        INFO("Primer name=" << read.name() << " " << ans);
    }
};

config::debruijn_config::pacbio_processor InitializePacBioProcessor() {
    config::debruijn_config::pacbio_processor pb;  
    pb.bwa_length_cutoff = 0; //500
    pb.compression_cutoff = 0.6; // 0.6
    pb.path_limit_stretching = 1.3; //1.3
    pb.path_limit_pressing = 0.7;//0.7
    pb.max_path_in_dijkstra = 15000; //15000
    pb.max_vertex_in_dijkstra = 15000; //2000
//gap_closer
    pb.long_seq_limit = 400; //400
    pb.pacbio_min_gap_quantity = 2; //2
    pb.contigs_min_gap_quantity = 1; //1
    pb.max_contigs_gap_length = 10000; // 10000
    return pb;
}

void LoadReads(const string &sequence_fasta, std::vector<io::SingleRead> &wrappedreads) {
    io::ReadStreamList<io::SingleRead> streams;
    streams.push_back(make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(sequence_fasta)));
    io::SingleStreamPtr sstream = io::MultifileWrap(streams); 
     
    while (!sstream->eof()) {
        io::SingleRead read;
        *sstream >> read;
        wrappedreads.push_back(std::move(read));
    }
    INFO("Loaded sequences from " << sequence_fasta);
}

void Launch(size_t K, const string &saves_path, const string &primer_fasta, const string &sequence_fasta, const string &mapper_type, const string &output_file, int threads) {
    conj_graph_pack gp(K, "tmp3", 0);
    graphio::ScanGraphPack(saves_path, gp);
    INFO("Loaded graph with " << gp.g.size() << " vertices");

    std::vector<io::SingleRead> wrappedprimers;
    std::vector<io::SingleRead> wrappedreads;
    LoadReads(primer_fasta, wrappedprimers);
    //LoadReads(sequence_fasta, wrappedreads);
    
    alignment::BWAIndex::AlignmentMode mode;
    if (mapper_type == "pacbio"){
        mode = alignment::BWAIndex::AlignmentMode::PacBio;
    } else if (mapper_type == "nanopore"){
        mode = alignment::BWAIndex::AlignmentMode::Ont2D;
    } else if (mapper_type == "16S"){
        mode = alignment::BWAIndex::AlignmentMode::Rna16S;
    } else {
        mode = alignment::BWAIndex::AlignmentMode::Default;
    }

    config::debruijn_config::pacbio_processor pb = InitializePacBioProcessor();
    PacBioAligner aligner(gp, mode, pb, output_file); 
    INFO("PacBioAligner created");

    ofstream myfile;
    myfile.open(output_file + ".gpa", std::ofstream::out | std::ofstream::app);
    myfile << "H\n";
    myfile.close();

    std::vector<ReadMapping> reads;
// #pragma omp parallel num_threads(threads)
// #pragma omp for
    for (size_t i =0 ; i < wrappedprimers.size(); ++i) {
        std::vector<int> v(50);
        std::vector<EdgeId> e(50);
        std::vector<MappingRange> range;
        aligner.AlignPrimer(wrappedprimers[i], v, e, range);
        reads.push_back(ReadMapping(wrappedprimers[i], v, e, range));
    }   

// #pragma omp parallel num_threads(threads)
// #pragma omp for
//         for (size_t i =0 ; i < wrappedreads.size(); ++i) {
//             aligner.AlignRead(wrappedreads[i]);
//         }
}
}

int main(int argc, char **argv) {
    if (argc < 7) {
        cout << "Usage: longreads_aligner <K>"
             << " <saves path> <primer file(fasta/fastq)> <sequences file (fasta/fastq)> <mapper type: {pacbio, nanopore, default} > <ouput-prefix> {threads-num(16)}" << endl;
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
    string primer_file = argv[3];
    INFO("Load primers from " << primer_file);
    string sequence_file = argv[4];
    INFO("Load sequences from " << sequence_file);
    string mapper_type = argv[5];
    INFO("Mapper type " << mapper_type);
    string output_file = argv[6];
    debruijn_graph::Launch(K, saves_path, primer_file, sequence_file, mapper_type, output_file, threads);
    return 0;
}
