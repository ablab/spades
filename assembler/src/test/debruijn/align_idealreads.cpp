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

public:
    IdealAligner(const conj_graph_pack &gp, 
                 const string &output_file,
                 std::shared_ptr<MapperClass> mapper):
      gp_(gp), output_file_(output_file), mapper_(mapper){}


    
    void AlignRead(const io::SingleRead &read){
        Sequence seq(read.sequence());
        INFO("Read " << read.name() <<". Current Read")
        auto current_mapping = mapper_->MapRead(read);
        if (current_mapping.empty()){
            INFO("Read " << read.name() <<" wasn't aligned");
        }
        if (current_mapping.size() > 0){
            std::string pathStr = "";
            std::string pathlenStr = "";
            for (const auto &edge_mapping: current_mapping) {
                EdgeId edgeid = edge_mapping.first;
                omnigraph::MappingRange range= edge_mapping.second;
                VertexId v1 = gp_.g.EdgeStart(edgeid);
                VertexId v2 = gp_.g.EdgeEnd(edgeid);
                //pathStr += std::to_string(edgeid.int_id()) + " (" + std::to_string(v1.int_id()) + "," + std::to_string(v2.int_id()) + ") ";
                pathStr += std::to_string(edgeid.int_id()) + ",";
                pathlenStr += std::to_string(range.mapped_range.end_pos - range.mapped_range.start_pos) + ",";
            }
            INFO("Path: " << pathStr);
            INFO("Read " << read.name() << " length=" << seq.size() << "; path_len=" << current_mapping.size()  << "; aligned: " << pathStr);
            //INFO("Seq subs: " << subStr);
#pragma omp critical
        {
            ofstream myfile;
            myfile.open(output_file_ + ".tsv", std::ofstream::out | std::ofstream::app);
            myfile << read.name() << "\t" << seq.size() << "\t" << current_mapping.size() << "\t" << pathStr << "\t" << pathlenStr << "\n";
            myfile.close();
        }
        }  
    }  
};


void Launch(size_t K, const string &saves_path, const string &reads_fasta, const string &output_file) {
    conj_graph_pack gp(K, "tmp3", 0);
    //graphio::ScanGraphPack(saves_path, gp);
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
        for (size_t i =0 ; i < wrappedreads.size(); ++i) {
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
