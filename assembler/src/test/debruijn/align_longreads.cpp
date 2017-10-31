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
      gp_(gp),pac_index_(gp.g, pb, mode), output_file_(output_file){}



    std::string print(map<string, string> &line) {
        std::vector<string> v = {"Ind", "Name", "ReadName", "StartR", "LenR", "DirR", "EdgeId", "StartE", "LenE", "DirE", "CIGAR", "Prev", "Next"};
        string outStr = "";
        for (const auto &it : v){
            outStr += line[it] + "\t";
        }
        return outStr;
    }


    void ToGPA(const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &aligned_mappings, const io::SingleRead &read) {
        int nameIndex = 0;
        std::string res = "";
        for (const auto &mappingpath : aligned_mappings){
            string prev = "";
            for (int i = 0; i < (int) mappingpath.size(); ++ i) {
                EdgeId edgeid = mappingpath.edge_at(i);
                omnigraph::MappingRange mapping = mappingpath.mapping_at(i);
                map<string, string> line = {{"Ind", "A"}, {"Name", ""}, {"ReadName", read.name()}, {"StartR", ""}, {"LenR", ""}, {"DirR", ""}
                                                                      , {"EdgeId", ""}, {"StartE", ""}, {"LenE", ""}, {"DirE", ""}
                                                                      , {"CIGAR", ""}, {"Prev", ""} , {"Next", ""}};
                
                nameIndex ++;                
                line["Name"] = read.name() + "_" + std::to_string(nameIndex);
                
                line["StartR"] = std::to_string(mapping.initial_range.start_pos); // TODO
                line["LenR"] = std::to_string(mapping.initial_range.end_pos - mapping.initial_range.start_pos); // TODO
                line["DirR"] = "?"; // TODO


                line["EdgeId"] = std::to_string(edgeid.int_id());
                line["StartE"] = std::to_string(mapping.mapped_range.start_pos);
                line["LenE"] = std::to_string(mapping.mapped_range.end_pos - mapping.mapped_range.start_pos);
                line["DirE"] = "?"; // TODO

                line["CIGAR"] = "*";// TODO mapping.cigar;

                if (i > 0){
                    line["Prev"] = prev;
                } else {
                    line["Prev"] = "-";
                }
                prev = line["Name"]; 
                if (i < (int) mappingpath.size() - 1){
                    line["Next"] = read.name() + "_" + std::to_string(nameIndex + 1);
                } else {
                    line["Next"] = "-";
                }
                res += print(line) + "\n";
                
            }
        }
#pragma omp critical
        {
            ofstream myfile;
            myfile.open(output_file_ + ".gpa", std::ofstream::out | std::ofstream::app);
            myfile << res;
            myfile.close();
        }        
    }
    
    void AlignRead(const io::SingleRead &read){
        Sequence seq(read.sequence());
        INFO("Read " << read.name() <<". Current Read")
        auto current_read_mapping = pac_index_.GetReadAlignment(seq);
        const auto& aligned_mappings = current_read_mapping.main_storage;
        if (aligned_mappings.size() > 0){
           // INFO("Read " << read.name() << " aligned with score=" << score << " and length=" << read.sequence().size());
        } else {
            INFO("Read " << read.name() << " wasn't aligned and length=" << read.sequence().size());
        }
        
        if (aligned_mappings.size() > 0){
            ToGPA(aligned_mappings, read);
            string pathStr = "";
            for (const auto &mappingpath : aligned_mappings){
                for (const auto &edgeid: mappingpath.simple_path()) {
                    VertexId v1 = gp_.g.EdgeStart(edgeid);
                    VertexId v2 = gp_.g.EdgeEnd(edgeid);
                    pathStr += std::to_string(edgeid.int_id()) + " (" + std::to_string(v1.int_id()) + "," + std::to_string(v2.int_id()) + ") ";
                }
                pathStr += "\n";
            }
            INFO("Paths: " << pathStr);
            pathStr = "";
            string subStr = "";
            string max_str = "";
            int max_len = 0;
            string sum_str = "";
            for (const auto &path : aligned_mappings){
                int seq_start = -1;
                int seq_end = 0;
                int mapping_start = 0;
                int mapping_end = 0;
                string cur_path = "";
                string cur_path_len = "";
                string cur_substr = "";
                string str = "";
                EdgeId last_edge= EdgeId();
                omnigraph::MappingRange last_range;
                for (int i = 0; i < (int) path.size(); ++ i) {
                    EdgeId edgeid = path.edge_at(i);
                    omnigraph::MappingRange mapping = path.mapping_at(i);
                    mapping_start = mapping.mapped_range.start_pos;
                    mapping_end = mapping.mapped_range.end_pos;
                    if (i > 0){
                        mapping_start = 0;
                    }
                    if (i < path.size() - 1) {
                        mapping_end = gp_.g.length(edgeid);
                    }
                    last_edge = edgeid;
                    last_range = mapping;
                    cur_path += std::to_string(edgeid.int_id()) + " (" + std::to_string(mapping_start) + "," + std::to_string(mapping_end) + ") ["
                               + std::to_string(mapping.initial_range.start_pos) + "," + std::to_string(mapping.initial_range.end_pos) + "], ";
                    
                    string tmp = gp_.g.EdgeNucls(edgeid).str();
                    //cur_path += std::to_string(mapping.edgeId.int_id()) + ",";
                    cur_path_len += std::to_string(mapping_end - mapping_start) + ",";
                    cur_substr += std::to_string(mapping.initial_range.start_pos) + "-" + std::to_string(mapping.initial_range.end_pos) + ", ";
                    string to_add = tmp.substr(mapping_start, mapping_end - mapping_start);
                    str += to_add;
                    if (seq_start < 0){
                        seq_start = (int) mapping.initial_range.start_pos;
                    }
                    seq_end = (int) mapping.initial_range.end_pos;
                }
                string tmp = gp_.g.EdgeNucls(last_edge).str();
                str += tmp.substr(last_range.mapped_range.end_pos, gp_.g.k());
                pathStr += cur_path + "; ";
                subStr += cur_substr + "\n";
                if (seq_end - seq_start > max_len){
                    max_len = seq_end - seq_start;
                    max_str = read.name() + "\t" + std::to_string(seq_start) + "\t" + std::to_string(seq_end) + "\t"  
                                                    + std::to_string(read.sequence().size())+  "\t" + cur_path + "\t" + cur_path_len + "\t" + str + "\n";
                }
                sum_str += read.name() + "\t" + std::to_string(seq_start) + "\t" + std::to_string(seq_end) + "\t" 
                                                 + std::to_string(read.sequence().size())+  "\t" + cur_path + "\t" + cur_path_len + "\t" + str + "\n";
            }
            INFO("Read " << read.name() << " aligned and length=" << read.sequence().size() <<  " and max_len=" << max_len);
            INFO("Read " << read.name() << ". Paths with ends: " << pathStr );
            //INFO("Seq subs: " << subStr);
#pragma omp critical
        {
            ofstream myfile;
            myfile.open(output_file_ + ".tsv", std::ofstream::out | std::ofstream::app);
            myfile << max_str;
            myfile.close();
        }
        }  
    }  
};

config::debruijn_config::pacbio_processor InitializePacBioProcessor() {
    config::debruijn_config::pacbio_processor pb;  
    pb.bwa_length_cutoff = 500; //500
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

void Launch(size_t K, const string &saves_path, const string &sequence_fasta, const string &mapper_type, const string &output_file) {
    conj_graph_pack gp(K, "tmp3", 0);
    graphio::ScanGraphPack(saves_path, gp);
    INFO("Loaded graph with " << gp.g.size() << " vertices");

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

#pragma omp parallel num_threads(16)
#pragma omp for
        for (size_t i =0 ; i < wrappedreads.size(); ++i) {
            aligner.AlignRead(wrappedreads[i]);
        }
}
}

int main(int argc, char **argv) {
    omp_set_num_threads(16);
    if (argc < 6) {
        cout << "Usage: longreads_aligner <K>"
             << " <saves path> <sequences file (fasta/fastq)> <mapper type: {pacbio, nanopore, default} > <ouput-prefix>" << endl;
        exit(1);
    }
    create_console_logger();
    size_t K = std::stoll(argv[1]);
    string saves_path = argv[2];
    INFO("Load graph from " << saves_path);
    string sequence_file = argv[3];
    INFO("Load sequences from " << sequence_file);
    string mapper_type = argv[4];
    INFO("Mapper type " << mapper_type);
    string output_file = argv[5];
    debruijn_graph::Launch(K, saves_path, sequence_file, mapper_type, output_file);
    return 0;
}
