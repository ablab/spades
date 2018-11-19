//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <iostream>
#include <fstream>

#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"

#include "io/binary/graph.hpp"
#include "io/reads/io_helper.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "io/reads/multifile_reader.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/graph/gfa_writer.hpp"

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "assembly_graph/stats/picture_dump.hpp"

using namespace std;

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

typedef debruijn_graph::BasicSequenceMapper<debruijn_graph::Graph, Index> MapperClass;

string getStrId(const EdgeId &e, const debruijn_graph::ConjugateDeBruijnGraph &g_) {
    if (e.int_id() < g_.conjugate(e).int_id()) {
        return to_string(e.int_id()) + "+";
    } else {
        return to_string(e.int_id()) + "-";
    }
}

class IdealAligner {
  private:
    const debruijn_graph::ConjugateDeBruijnGraph &g_;
    const string &output_file_;
    shared_ptr<MapperClass> mapper_;
    ofstream myfile_;

  public:
    IdealAligner(const debruijn_graph::ConjugateDeBruijnGraph &g,
                 const string &output_file,
                 shared_ptr<MapperClass> mapper):
        g_(g), output_file_(output_file), mapper_(mapper) {
        myfile_.open(output_file_ + ".tsv", ofstream::out);
    }


    bool IsCanonical(EdgeId e) const {
        return e <= g_.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : g_.conjugate(e);
    }

    vector<int> CheckPathConsistency(const vector<EdgeId> &read_mapping, const MappingPath<EdgeId>& current_mapping)  {
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

    void AlignRead(const io::SingleRead &read) {
        Sequence seq(read.sequence());
        INFO("Read " << read.name() << ".")
        auto current_mapping = mapper_->MapRead(read);
        ReadPathFinder<debruijn_graph::Graph> readmapper(g_);
        auto read_mapping = readmapper.FindReadPath(current_mapping);
        if (current_mapping.empty() || read_mapping.size() == 0) {
            INFO("Read " << read.name() << " wasn't aligned");
        }

        vector<int> inds = CheckPathConsistency(read_mapping, current_mapping);
        if (inds.size() > 0) {
            int j = 0;
            int len_before = 0;
            string path_str = "";
            string pathlen_str = "";
            string edgelen_str = "";
            string str = "";
            for (int i = 0; i < (int) read_mapping.size(); ++ i) {
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
                    if ( (i > 0 && i < (int)(read_mapping.size()) - 1) && (mapping_end - mapping_start != initial_end - initial_start || mapping_end - mapping_start != g_.length(edgeid)) ) {
                        DEBUG("Bad ranges")
                        return;
                    }
                    ++ j;
                }
                len_before += (int) g_.length(edgeid);
                path_str += getStrId(edgeid, g_) + " (" + to_string(mapping_start) + "," + to_string(mapping_end) + ") ["
                            + to_string(initial_start) + "," + to_string(initial_end) + "], ";
                pathlen_str += to_string(mapping_end - mapping_start) + ",";
                edgelen_str += to_string(g_.length(edgeid)) + ",";
                str += g_.EdgeNucls(edgeid).Subseq(mapping_start, mapping_end).str();
            }
            DEBUG("Path: " << path_str);
            DEBUG("Read " << read.name() << " length=" << seq.size() << "; path_len=" << current_mapping.size()  << "; aligned: " << path_str);
            if (str.size() != seq.size()) {
                DEBUG("Read " << read.name() << " wasn't fully aligned");
                return;
            }
            #pragma omp critical
            {
                myfile_ << read.name() << "\t" << seq.size() << "\t" << current_mapping.size()
                        << "\t" << path_str << "\t" << pathlen_str << "\t" << edgelen_str << "\t" << str << "\n";

            }
        }
    }

};

void LoadGraph(const string &saves_path, bool load_spades_graph, debruijn_graph::ConjugateDeBruijnGraph &g) {
    if (fs::extension(saves_path) == ".gfa") {
        DEBUG("Load gfa")
        VERIFY_MSG(fs::is_regular_file(saves_path), "GFA-file " + saves_path + " doesn't exist");
        gfa::GFAReader gfa(saves_path);
        DEBUG("Segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
        gfa.to_graph(g, load_spades_graph);
        return;
    } else {
        INFO("Load from saves");
        io::binary::Load(saves_path, g);
        return;
    }
}


void Launch(size_t K, const string &saves_path, bool load_spades_graph, const string &reads_fasta, const string &output_file) {
    string tmpdir = fs::make_temp_dir(fs::current_dir(), "tmp");
    debruijn_graph::ConjugateDeBruijnGraph g(K);
    LoadGraph(saves_path, load_spades_graph, g);

    KmerMapper<debruijn_graph::Graph> kmer_mapper(g);
    kmer_mapper.Detach();
    EdgeIndex<debruijn_graph::Graph> edge_index(g, tmpdir);
    edge_index.Detach();
    shared_ptr<MapperClass> mapper(new MapperClass(g, edge_index, kmer_mapper));
    edge_index.Attach();
    kmer_mapper.Attach();
    edge_index.Refill();
    INFO("Loaded graph with " << g.size() << " vertices");

    io::ReadStreamList<io::SingleRead> streams;
    streams.push_back(make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(reads_fasta)));

    io::SingleStreamPtr sstream = io::MultifileWrap(streams);
    vector<io::SingleRead> wrappedreads;
    while (!sstream->eof()) {
        io::SingleRead read;
        *sstream >> read;
        wrappedreads.push_back(move(read));
    }
    INFO("Loaded reads from " << reads_fasta);

    IdealAligner aligner(g, output_file, mapper);
    INFO("IdealAligner created");

    #pragma omp parallel num_threads(16)
    #pragma omp for
    for (size_t i = 0 ; i < wrappedreads.size(); ++i) {
        aligner.AlignRead(wrappedreads[i]);
    }
    if (!load_spades_graph) {
        INFO("Saving *.gfa to " << output_file + ".gfa");
        ofstream os(output_file + ".gfa");
        gfa::GFAWriter gfa_writer(g, os);
        gfa_writer.WriteSegmentsAndLinks();
    }
    fs::remove_dir(tmpdir);
}
}

int main(int argc, char **argv) {
    omp_set_num_threads(1);
    if (argc < 5) {
        cout << "Usage: idealreads_aligner <K>"
             << " <saves path> <long reads file (fasta)> <ouput-prefix> --spades" << endl;
        exit(1);
    }
    create_console_logger();
    size_t K = stoll(argv[1]);
    string saves_path = argv[2];
    INFO("Load graph from " << saves_path);
    string longreads_file = argv[3];
    INFO("Load long reads from " << longreads_file);
    string output_file = argv[4];
    bool load_spades_graph = false;
    if (argc == 6) {
        load_spades_graph = true;
    }
    debruijn_graph::Launch(K, saves_path, load_spades_graph, longreads_file, output_file);
    return 0;
}
