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

//#include "modules/alignment/pacbio/gap_dijkstra.hpp"
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

    ReadMapping(const io::SingleRead &read, int &v, EdgeId &e, MappingRange &range) 
        : read_(read), v_(v), e_(e), range_(range) {
        ExtractPositions();
    }

    bool operator < (const ReadMapping &mapping) const {
        return (this->start_pos_ < mapping.start_pos_ || (this->start_pos_ == mapping.start_pos_ && this->end_pos_ < mapping.start_pos_) );
    }

    io::SingleRead read_;
    int v_;
    EdgeId e_;
    MappingRange range_;
    int start_pos_;
    int end_pos_;
};


class SequenceAligner {
private:    
    const conj_graph_pack &gp_;
    const pacbio::PacBioMappingIndex<Graph> pac_index_;
    const string &output_file_;

    std::vector<ReadMapping> primers_;

public:
    SequenceAligner(const conj_graph_pack &gp, 
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

    void AlignPrimer(const io::SingleRead &read, int &v, EdgeId &e, MappingRange &range) {
        v = -1;
        e = EdgeId();
        range = MappingRange(-1, -1, -1, -1);
        for (auto it = gp_.g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId eid = *it;
            std::string edge_str = gp_.g.EdgeNucls(eid).str();
            int start_pos = -1;
            int end_pos = -1;
            int dist = EditDistance(read.sequence().str(), edge_str, start_pos, end_pos);
            if (dist == -1) continue;
            if (v == -1 || v > dist) {
                v = dist;
                e = eid;
                range = MappingRange(Range(0, read.size()), Range(start_pos, end_pos + 1));
            }
        }
        std::string ans = " best: edge_id=" + std::to_string(e.int_id()) + " dist=" +  std::to_string(v);
        INFO("Primer name=" << read.name() << " " << ans);
    }

    void PreparePrimers(std::vector<io::SingleRead> &wrappedprimers) {
        // #pragma omp parallel num_threads(threads)
        // #pragma omp for
        for (size_t i =0 ; i < wrappedprimers.size(); ++i) {
            int v;
            EdgeId e;
            MappingRange range;
            AlignPrimer(wrappedprimers[i], v, e, range);
            primers_.push_back(ReadMapping(wrappedprimers[i], v, e, range));
        }   
        std::sort(primers_.begin(), primers_.end());
    }

    void PrepareInitialState(omnigraph::MappingPath<debruijn_graph::EdgeId> &path, const Sequence &s, bool forward, Sequence &ss, EdgeId &start_e, int &start_pos) const {
        if (forward){
            start_e = path.edge_at(path.size() - 1);
            omnigraph::MappingRange mapping = path.mapping_at(path.size() - 1);
            start_pos = mapping.mapped_range.end_pos;
            ss = s.Subseq(mapping.initial_range.end_pos, (int) s.size() );
        } else {
            start_e = gp_.g.conjugate(path.edge_at(0));
            omnigraph::MappingRange mapping = path.mapping_at(0);
            start_pos = min(gp_.g.length(start_e), gp_.g.length(start_e) + gp_.g.k() - mapping.mapped_range.start_pos);
            ss = !s.Subseq(0, mapping.initial_range.start_pos);
        }
    }

    void UpdatePath(omnigraph::MappingPath<debruijn_graph::EdgeId> &path, std::vector<EdgeId> &ans, int end_pos, bool forward) const {
        if (forward) {
            for (int i = 1; i < ans.size() - 1; ++i) {
                path.push_back(ans[i], omnigraph::MappingRange(Range(0, 0), Range(0, gp_.g.length(ans[i])) ));
            }
            path.push_back(ans[ans.size() - 1], omnigraph::MappingRange(Range(0, 0), Range(0, end_pos -  gp_.g.k()) ));
        } else {
            omnigraph::MappingPath<debruijn_graph::EdgeId> cur_sorted;
            int start = gp_.g.length(ans[ans.size() - 1]) + gp_.g.k() - end_pos;
            int cur_ind = ans.size() - 1;
            while (cur_ind >= 0 && start - (int) gp_.g.length(ans[cur_ind]) > 0){
                start -= gp_.g.length(ans[cur_ind]);
                cur_ind --;
            }
            if (cur_ind > 0){
                cur_sorted.push_back(gp_.g.conjugate(ans[cur_ind]), omnigraph::MappingRange(Range(0, 0), Range(start, gp_.g.length(ans[cur_ind])) ));
            }
            for (int i = cur_ind - 1; i > 0; --i) {
                cur_sorted.push_back(gp_.g.conjugate(ans[i]), omnigraph::MappingRange(Range(0, 0), Range(0, gp_.g.length(ans[i])) ));
            }
            for (int i = 0; i < path.size(); ++i) {
                cur_sorted.push_back(path[i].first, path[i].second);
            }
            path = cur_sorted;
        }
    }

    void GrowEnds(omnigraph::MappingPath<debruijn_graph::EdgeId> &path, const Sequence &s, bool forward) const {
        VERIFY(path.size() > 0);
        Sequence ss; 
        int start_pos = -1;
        EdgeId start_e = EdgeId();
        PrepareInitialState(path, s, forward, ss, start_e, start_pos);

        int s_len = int(ss.size());
        int score = max(20, s_len/4);
        if (s_len > 2000) {
            INFO("EdgeDijkstra: sequence is too long " << s_len)
            return;
        }
        if (s_len < gp_.g.length(start_e) + gp_.g.k() - start_pos) {
            INFO("EdgeDijkstra: sequence is too small " << s_len)
            return;
        }
        pacbio::DijkstraEndsReconstructor algo = pacbio::DijkstraEndsReconstructor(gp_.g, ss, start_e, start_pos, score);
        algo.CloseGap();
        score = algo.GetEditDistance();
        if (score == -1){
            INFO("EdgeDijkstra didn't find anything edge=" << start_e.int_id() << " s_start=" << start_pos << " seq_len=" << ss.size())
            return;
        }

        std::vector<EdgeId> ans = algo.GetPath();
        int end_pos = algo.GetPathEndPosition();
        UpdatePath(path, ans, end_pos, forward);
    }

    std::vector<EdgeId> FillGap(MappingRange &a_range, EdgeId &a_e, MappingRange &b_range, EdgeId &b_e, const Sequence &s) {
        auto path_searcher_b = omnigraph::DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(g_, path_max_length);
        path_searcher_b.Run(end_v);
        auto path_searcher = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, path_max_length);
        path_searcher.Run(start_v);
        auto reached_vertices_b = path_searcher_b.ProcessedVertices();
        auto reached_vertices = path_searcher.ProcessedVertices();

        std::map<VertexId, size_t> vertex_pathlen;
        for (auto j_iter = reached_vertices_b.begin(); j_iter != reached_vertices_b.end(); ++j_iter) {
                if (reached_vertices.count(*j_iter) > 0){
                        vertex_pathlen[*j_iter] = path_searcher_b.GetDistance(*j_iter);
                }
        }
        if (end_pos < start_pos) {
            WARN ("modifying limits because of some bullshit magic, seq length 0")
            end_pos = start_pos;
        }

        Sequence ss = s.Subseq(start_pos, min(end_pos + 1, int(s.size()) ));
        int s_len = int(ss.size());
        path_max_length = min(score, max(s_len/3, 20));
        DEBUG(" Dijkstra: String length " << s_len << " max-len " << path_max_length << " start_p=" << start_p << " end_p=" << end_p);
        if (s_len > 2000 && vertex_pathlen.size() > 100000){
            INFO("Dijkstra won't run: Too big gap or too many paths");
            return vector<EdgeId>(0);
        }
        DijkstraGapFiller gap_filler = DijkstraGapFiller(g_, ss, start_e, end_e, start_p, end_p, path_max_length, vertex_pathlen);
        gap_filler.CloseGap();
        score = gap_filler.GetEditDistance();
        if (score == -1){
            INFO("Dijkstra didn't find anything")
            score = STRING_DIST_INF;
            return vector<EdgeId>(0);
        }
        std::vector<EdgeId> ans = gap_filler.GetPath();        

        return ans;
    }

    omnigraph::MappingPath<debruijn_graph::EdgeId> FillGapsBetweenPrimers(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappings, const io::SingleRead &read) {
        omnigraph::MappingPath<debruijn_graph::EdgeId> res;
        res.push_back(mappings.edge_at(0), mappings.mapping_at(0));
        GrowEnds(res, read.sequence(), false);
        for (size_t i = 0; i < mappings.size(); ++ i){
            if (i + 1 < mappings.size()) {
                if (Consistent(mappings[i], mappings[i + 1])) {
                    if (Adjacent(mappings[i], mappings[i + 1])) {
                        res.push_back(mappings.edge_at(i + 1), mappings.mapping_at(i + 1));
                    } else {
                        std::vector<EdgeId> path = FillGap(mappings[i], mappings[i + 1], read.sequence());   
                    }
                }
            }
        }

    }

    void AlignRead(const io::SingleRead &read) {
        omnigraph::MappingPath<debruijn_graph::EdgeId> mappings;
        int threshold = 20;
        for (const auto primer: primers_) {
            int start_pos = -1;
            int end_pos = -1;
            int dist = EditDistance(primer.read_.sequence().str(), read.sequence().str(), start_pos, end_pos);
            if (dist != -1 && abs(start_pos - primer.start_pos_) < threshold && abs(end_pos - primer.end_pos_) < threshold){
                mappings.push_back(primer.e_, MappingRange(Range(start_pos, end_pos), 
                                                          Range(primer.range_.mapped_range.start_pos, primer.range_.mapped_range.start_pos)));
            }
        }
        if (mappings.size() > 0) {
            omnigraph::MappingPath<debruijn_graph::EdgeId> res = FillGapsBetweenPrimers(mappings, read);
        } else {
            omnigraph::MappingPath<debruijn_graph::EdgeId>();
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
    SequenceAligner aligner(gp, mode, pb, output_file); 
    INFO("SequenceAligner created");

    ofstream myfile;
    myfile.open(output_file + ".gpa", std::ofstream::out | std::ofstream::app);
    myfile << "H\n";
    myfile.close();

    aligner.PreparePrimers(wrappedprimers);

#pragma omp parallel num_threads(threads)
#pragma omp for
    for (size_t i =0 ; i < wrappedreads.size(); ++i) {
        aligner.AlignRead(wrappedreads[i]);
    }
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
