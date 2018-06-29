//***************************************************************************
//* Copyright (c) 2015-2018 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/stl_utils.hpp"
#include "utils/logger/log_writers.hpp"

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
#include "pipeline/config_struct.hpp"

#include "assembly_graph/core/graph.hpp"
#include "modules/alignment/sequence_mapper.hpp"

#include "io/reads/io_helper.hpp"
#include "common/assembly_graph/core/coverage.hpp"

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

typedef debruijn_graph::BasicSequenceMapper<debruijn_graph::Graph, Index> MapperClass;

class IdealAligner {
  private:
    const conj_graph_pack &gp_;
    const string &output_file_;
    std::shared_ptr<MapperClass> mapper_;
    ofstream myfile_;

  public:
    IdealAligner(const conj_graph_pack &gp,
                 const string &output_file,
                 std::shared_ptr<MapperClass> mapper):
        gp_(gp), output_file_(output_file), mapper_(mapper) {
        myfile_.open(output_file_ + ".tsv", std::ofstream::out);
    }


    bool IsCanonical(EdgeId e) const {
        return e <= gp_.g.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : gp_.g.conjugate(e);
    }

    std::vector<int> CheckPathConsistency(const vector<EdgeId> &read_mapping, const MappingPath<EdgeId>& current_mapping)  {
        std::vector<int> inds;
        int j = 0;
        int len_before = 0;
        for (int i = 0; i < read_mapping.size(); ++ i) {
            if (i > 0) {
                VertexId v1 = gp_.g.EdgeEnd(read_mapping[i - 1]);
                VertexId v2 = gp_.g.EdgeStart(read_mapping[i]);
                if (v1 != v2) {
                    INFO("Not a connected path!")
                    return std::vector<int>();
                }
            }
            EdgeId edgeid = read_mapping[i];
            if (j < current_mapping.size() && current_mapping[j].first == edgeid) {
                omnigraph::MappingRange range = current_mapping[j].second;
                if (len_before == range.initial_range.start_pos) {
                    inds.push_back(i);
                    j++;
                    len_before = range.initial_range.end_pos;
                } else {
                    len_before += gp_.g.length(edgeid);    
                }
            } else {
                len_before += gp_.g.length(edgeid);
            }

        }
        if (j != current_mapping.size()) {
            INFO("Badly mapped j != mapping")
            return std::vector<int>();
        }
        return inds;
    }

    void AlignRead(const io::SingleRead &read) {
        Sequence seq(read.sequence());
        INFO("Read " << read.name() << ". Current Read")
        auto current_mapping = mapper_->MapRead(read);
        ReadPathFinder<debruijn_graph::Graph> readmapper(gp_.g);
        auto read_mapping = readmapper.FindReadPath(current_mapping);
        if (current_mapping.empty() || read_mapping.size() == 0) {
            INFO("Read " << read.name() << " wasn't aligned");
        }

        std::vector<int> inds = CheckPathConsistency(read_mapping, current_mapping);
        if (inds.size() > 0){
            int j = 0;
            int len_before = 0;
            std::string path_str = "";
            std::string pathlen_str = "";
            std::string edgelen_str = "";
            std::string str = "";
            for (int i = 0; i < read_mapping.size(); ++ i) {
                EdgeId edgeid = read_mapping[i];
                size_t mapping_start = 0;
                size_t mapping_end = gp_.g.length(edgeid);
                size_t initial_start = len_before;
                size_t initial_end = len_before + gp_.g.length(edgeid);
                if (inds[j] == i) {
                    omnigraph::MappingRange range = current_mapping[j].second;
                    mapping_start = range.mapped_range.start_pos;
                    mapping_end = j + 1 < inds.size() ? range.mapped_range.end_pos : range.mapped_range.end_pos + gp_.g.k();
                    initial_start = range.initial_range.start_pos;
                    initial_end = j + 1 < inds.size() ? range.initial_range.end_pos : range.initial_range.end_pos + gp_.g.k();
                    if ( (i > 0 && i < read_mapping.size() - 1) && (mapping_end - mapping_start != initial_end - initial_start || mapping_end - mapping_start != gp_.g.length(edgeid)) ) {
                        INFO("Bad ranges")
                        return;
                    }
                    ++ j;
                }
                len_before += gp_.g.length(edgeid);
                path_str += std::to_string(edgeid.int_id()) + " (" + std::to_string(mapping_start) + "," + std::to_string(mapping_end) + ") ["
                            + std::to_string(initial_start) + "," + std::to_string(initial_end) + "], ";
                pathlen_str += std::to_string(mapping_end - mapping_start) + ",";
                edgelen_str += std::to_string(gp_.g.length(edgeid)) + ",";
                str += gp_.g.EdgeNucls(edgeid).Subseq(mapping_start, mapping_end).str();
            }
            INFO("Path: " << path_str);
            INFO("Read " << read.name() << " length=" << seq.size() << "; path_len=" << current_mapping.size()  << "; aligned: " << path_str);
            #pragma omp critical
            {
                myfile_ << read.name() << "\t" << seq.size() << "\t" << current_mapping.size()
                        << "\t" << path_str << "\t" << pathlen_str << "\t" << edgelen_str << "\t" << str << "\n";

            }
        }
    }

};


void Launch(size_t K, const string &saves_path, const string &reads_fasta, const string &output_file) {
    conj_graph_pack gp(K, "tmp3", 0);
    std::shared_ptr<MapperClass> mapper(new MapperClass(gp.g, gp.index, gp.kmer_mapper));
    gp.kmer_mapper.Attach();
    debruijn_graph::graphio::ScanGraphPack(saves_path, gp);
    INFO("Loaded graph with " << gp.g.size() << " vertices");

    io::ReadStreamList<io::SingleRead> streams;
    streams.push_back(make_shared<io::FixingWrapper>(make_shared<io::FileReadStream>(reads_fasta)));

    io::SingleStreamPtr sstream = io::MultifileWrap(streams);
    std::vector<io::SingleRead> wrappedreads;
    while (!sstream->eof()) {
        io::SingleRead read;
        *sstream >> read;
        wrappedreads.push_back(std::move(read));
    }
    INFO("Loaded reads from " << reads_fasta);

    IdealAligner aligner(gp, output_file, mapper);
    INFO("IdealAligner created");

    #pragma omp parallel num_threads(16)
    #pragma omp for
    for (size_t i = 0 ; i < wrappedreads.size(); ++i) {
        aligner.AlignRead(wrappedreads[i]);
    }
}
}

int main(int argc, char **argv) {
    omp_set_num_threads(16);
    if (argc < 5) {
        cout << "Usage: idealreads_aligner <K>"
             << " <saves path> <long reads file (fasta)> <ouput-prefix>" << endl;
        exit(1);
    }
    create_console_logger();
    size_t K = std::stoll(argv[1]);
    string saves_path = argv[2];
    INFO("Load graph from " << saves_path);
    string longreads_file = argv[3];
    INFO("Load long reads from " << longreads_file);
    string output_file = argv[4];
    debruijn_graph::Launch(K, saves_path, longreads_file, output_file);
    return 0;
}
