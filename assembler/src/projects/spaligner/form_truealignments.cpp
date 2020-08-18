//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"

#include "io/utils/edge_namer.hpp"
#include "io/binary/graph.hpp"
#include "io/reads/io_helper.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "io/reads/multifile_reader.hpp"
#include "io/reads/file_reader.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/graph/gfa_writer.hpp"

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "assembly_graph/stats/picture_dump.hpp"

#include <iostream>
#include <fstream>

using namespace std;

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

typedef BasicSequenceMapper<Graph, EdgeIndex<Graph>> MapperClass;


struct ReadMappingStr {
    string path_str = "";
    string pathlen_str = "";
    string edgelen_str = "";
    string str = "";
    size_t hits_num = 0;
    size_t paths_num = 0;
    bool mapped = false;

    ReadMappingStr(size_t hits_num_, size_t paths_num_)
    : hits_num(hits_num_),
      paths_num(paths_num_) {}
};

class KMerAligner {
  public:
    KMerAligner(const debruijn_graph::ConjugateDeBruijnGraph &g,
                 const io::CanonicalEdgeHelper<debruijn_graph::Graph> &edge_namer,
                 const string &output_file,
                 shared_ptr<MapperClass> mapper):
        g_(g), edge_namer_(edge_namer), output_file_(output_file), mapper_(mapper) {
        myfile_.open(output_file_ + ".tsv", ofstream::out);
    }

    ReadMappingStr AlignRead(const io::SingleRead &read) const {
        Sequence seq(read.sequence());
        INFO("Read " << read.name() << ".")
        auto current_mapping = mapper_->MapRead(read);
        ReadPathFinder<debruijn_graph::Graph> readmapper(g_);
        auto read_mapping = readmapper.FindReadPath(current_mapping);
        if (current_mapping.empty() || read_mapping.size() == 0) {
            INFO("Read " << read.name() << " wasn't aligned");
        }

        vector<int> inds = CheckPathConsistency(read_mapping, current_mapping);
        ReadMappingStr read_mapping_str(read_mapping.size(), current_mapping.size());
        if (inds.size() > 0) {
            int j = 0;
            int len_before = 0;
            for (int i = 0; i < (int) read_mapping.size(); ++ i) {
                if (!FormMappingString(read_mapping, inds, current_mapping, i, j, len_before, read_mapping_str)) {
                    return read_mapping_str;
                }
            }
            DEBUG("Path: " << read_mapping_str.path_str);
            DEBUG("Read " << read.name() << " length=" << seq.size() << "; path_len=" << current_mapping.size()  << "; aligned: " << read_mapping_str.path_str);
            if (read_mapping_str.str.size() != seq.size()) {
                DEBUG("Read " << read.name() << " wasn't fully aligned");
                return read_mapping_str;
            }
            read_mapping_str.mapped = true;
            return read_mapping_str;
        }
        return read_mapping_str;
    }

    void Print(const io::SingleRead &read, const ReadMappingStr &read_mapping_str){
        Sequence seq(read.sequence());
        #pragma omp critical
        {
            myfile_ << read.name() << "\t" << seq.size() << "\t" << read_mapping_str.hits_num
                    << "\t" << read_mapping_str.path_str << "\t" << read_mapping_str.pathlen_str 
                    << "\t" << read_mapping_str.edgelen_str << "\t" << read_mapping_str.str << "\n";

        }
    }

  private:

    string StrId(const EdgeId &e) const {
        return edge_namer_.EdgeOrientationString(e);
    }
    
    bool IsCanonical(EdgeId e) const {
        return e <= g_.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : g_.conjugate(e);
    }

    vector<int> CheckPathConsistency(const vector<EdgeId> &read_mapping, const MappingPath<EdgeId>& current_mapping) const {
        vector<int> inds;
        int j = 0;
        int len_before = 0;
        for (int i = 0; i < (int) read_mapping.size(); ++ i) {
            if (i > 0) {
                VertexId v1 = g_.EdgeEnd(read_mapping[i - 1]);
                VertexId v2 = g_.EdgeStart(read_mapping[i]);
                if (v1 != v2) {
                    DEBUG("Not a connected path!")
                    return vector<int>();
                }
            }
            EdgeId edgeid = read_mapping[i];
            if (j < (int) current_mapping.size() && current_mapping[j].first == edgeid) {
                omnigraph::MappingRange range = current_mapping[j].second;
                if (len_before == (int) range.initial_range.start_pos) {
                    inds.push_back(i);
                    j++;
                    len_before = (int) range.initial_range.end_pos;
                } else {
                    DEBUG("len_before: " << len_before << " start=" << range.initial_range.start_pos
                          << " end=" << range.initial_range.end_pos << " e_sz=" << g_.length(edgeid))
                    len_before += (int) g_.length(edgeid);
                }
            } else {
                len_before += (int) g_.length(edgeid);
            }

        }
        if (j != (int) current_mapping.size()) {
            DEBUG("Badly mapped j != mapping")
            return vector<int>();
        }
        return inds;
    }

    bool FormMappingString(const vector<EdgeId> &read_mapping, 
                           const vector<int> &inds,
                           const MappingPath<EdgeId>& current_mapping,
                           int i, int &j,
                           int &len_before, ReadMappingStr &read_mapping_str) const {
        EdgeId edgeid = read_mapping[i];
        size_t mapping_start = 0;
        size_t mapping_end = g_.length(edgeid);
        size_t initial_start = len_before;
        size_t initial_end = len_before + g_.length(edgeid);
        if (inds[j] == i) {
            omnigraph::MappingRange range = current_mapping[j].second;
            mapping_start = range.mapped_range.start_pos;
            mapping_end = j + 1 < (int) inds.size() ? range.mapped_range.end_pos : range.mapped_range.end_pos + g_.k();
            initial_start = range.initial_range.start_pos;
            initial_end = j + 1 < (int) inds.size() ? range.initial_range.end_pos : range.initial_range.end_pos + g_.k();
            if ((i > 0 && i < (int)(read_mapping.size()) - 1) &&
                    (mapping_end - mapping_start != initial_end - initial_start || mapping_end - mapping_start != g_.length(edgeid))) {
                DEBUG("Bad ranges")
                return false;
            }
            ++ j;
        }
        len_before += (int) g_.length(edgeid);
        read_mapping_str.path_str += StrId(edgeid) + " (" + to_string(mapping_start) + "," + to_string(mapping_end) + ") ["
                    + to_string(initial_start) + "," + to_string(initial_end) + "], ";
        read_mapping_str.pathlen_str += to_string(mapping_end - mapping_start) + ",";
        read_mapping_str.edgelen_str += to_string(g_.length(edgeid)) + ",";
        read_mapping_str.str += g_.EdgeNucls(edgeid).Subseq(mapping_start, mapping_end).str();
        return true;
    }

    const debruijn_graph::ConjugateDeBruijnGraph &g_;
    const io::CanonicalEdgeHelper<debruijn_graph::Graph> &edge_namer_;
    const string &output_file_;
    shared_ptr<MapperClass> mapper_;
    ofstream myfile_;


};

void LoadGraph(const string &saves_path, debruijn_graph::ConjugateDeBruijnGraph &g, io::IdMapper<std::string> &id_mapper) {
    if (fs::extension(saves_path) == ".gfa") {
        DEBUG("Load gfa")
        CHECK_FATAL_ERROR(fs::is_regular_file(saves_path), "GFA-file " + saves_path + " doesn't exist");
        gfa::GFAReader gfa(saves_path);
        DEBUG("Segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
        gfa.to_graph(g, &id_mapper);
        return;
    } else {
        INFO("Load from saves");
        io::binary::Load(saves_path, g);
        return;
    }
}


void Launch(size_t K, const string &saves_path, const string &reads_fasta, const string &output_file) {
    string tmpdir = fs::make_temp_dir(fs::current_dir(), "tmp");
    debruijn_graph::ConjugateDeBruijnGraph g(K);
    io::IdMapper<std::string> id_mapper;
    LoadGraph(saves_path, g, id_mapper);
    KmerMapper<debruijn_graph::Graph> kmer_mapper(g);
    kmer_mapper.Detach();
    EdgeIndex<debruijn_graph::Graph> edge_index(g, tmpdir);
    edge_index.Detach();
    shared_ptr<MapperClass> mapper(new MapperClass(g, edge_index, kmer_mapper));
    edge_index.Attach();
    kmer_mapper.Attach();
    edge_index.Refill();
    INFO("Loaded graph with " << g.size() << " vertices");

    auto sstream = io::FixingWrapper(io::FileReadStream(reads_fasta));
    std::vector<io::SingleRead> wrappedreads;
    while (!sstream.eof()) {
        io::SingleRead read;
        sstream >> read;
        wrappedreads.push_back(std::move(read));
    }
    INFO("Loaded reads from " << reads_fasta);

    io::CanonicalEdgeHelper<debruijn_graph::Graph> edge_namer(g, io::MapNamingF<debruijn_graph::Graph>(id_mapper));
    KMerAligner aligner(g, edge_namer, output_file, mapper);
    INFO("KMerAligner created");

    #pragma omp parallel num_threads(16)
    #pragma omp for
    for (size_t i = 0 ; i < wrappedreads.size(); ++i) {
        debruijn_graph::ReadMappingStr read_mapping_res = aligner.AlignRead(wrappedreads[i]);
        if (read_mapping_res.mapped) {
            aligner.Print(wrappedreads[i], read_mapping_res);
        }
    }
    fs::remove_dir(tmpdir);
}
}

int main(int argc, char **argv) {
    omp_set_num_threads(1);
    if (argc < 5) {
        cout << "Usage: idealreads_aligner <K>"
             << " <saves path> <long reads file (fasta)> <ouput-prefix>" << endl;
        exit(1);
    }
    create_console_logger();
    size_t K = stoll(argv[1]);
    string saves_path = argv[2];
    INFO("Load graph from " << saves_path);
    string longreads_file = argv[3];
    INFO("Load long reads from " << longreads_file);
    string output_file = argv[4];
    debruijn_graph::Launch(K, saves_path, longreads_file, output_file);
    return 0;
}
