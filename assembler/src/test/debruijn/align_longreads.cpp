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

    void getCIGAR(std::string &read, std::string aligned, std::string &cigar, int &score) {
        int d = max((int) read.size(), 20);
        edlib::EdlibAlignResult result = edlib::edlibAlign(aligned.c_str(), (int) aligned.size(), read.c_str(), (int) read.size()
                                           , edlib::edlibNewAlignConfig(d, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_PATH,
                                                                 NULL, 0));
        cigar = "";
        score = pacbio::STRING_DIST_INF;
        if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
            score = result.editDistance;
            cigar = edlib::edlibAlignmentToCigar(result.alignment, result.alignmentLength, edlib::EDLIB_CIGAR_EXTENDED);
        }
        edlib::edlibFreeAlignResult(result);
        string cur_num = "";
        int n = -1;
        int len_r = 0;
        int len_a = 0;
        for (size_t i = 0; i < cigar.size(); ++ i) {
            if (isdigit(cigar[i])) {
                cur_num += cigar[i];
            } else {
                n = std::stoi(cur_num);
                char c = cigar[i];
                if (c == '=' || c == 'I' || c == 'X' || c == 'M'){
                    len_a += n;
                }
                if (c != 'I'){
                    len_r += n;
                }
                cur_num = "";
                n = 0;
            }
        }
        DEBUG("CIGAR: "<< len_a << " " << aligned.size()  << " " << len_r << " " << read.size());
    }

    void getByEdgeCIGAR(string &read, string &aligned, std::vector<size_t> &edgeblocks, size_t start, std::vector<string> &edgecigar, std::vector<Range> &edge_initial_ranges, int &score) {
        std::string cigar;
        getCIGAR(read, aligned, cigar, score);
        DEBUG("CIGAR: " << "\n" << read << "\n" << aligned << "\n" << cigar );
        string cur_num = "";
        int n = 0;
        size_t r_i = 0;
        size_t a_i = 0;
        size_t cur_block = 0;
        string cur_cigar = "";
        size_t cur_start_pos = start;
        size_t cur_end_pos = start;
        for (size_t i = 0; i < cigar.size(); ++ i) {
            if (isdigit(cigar[i])) {
                cur_num += cigar[i];
            } else {
                n = std::stoi(cur_num);
                char c = cigar[i];
                if (c == '=' || c == 'I' || c == 'X' || c == 'M'){
                    while (a_i + n > edgeblocks[cur_block]) {
                        DEBUG("CIGAR: " << n << c);
                        n -= (int) (edgeblocks[cur_block] - a_i);
                        if (c != 'I') {
                            r_i += (edgeblocks[cur_block] - a_i);
                            cur_end_pos += (edgeblocks[cur_block] - a_i);
                        }
                        edge_initial_ranges.push_back(Range(cur_start_pos, cur_end_pos));
                        cur_start_pos = cur_end_pos;
                        if (edgeblocks[cur_block] - a_i != 0) {
                            edgecigar.push_back(cur_cigar + std::to_string(edgeblocks[cur_block] - a_i) + c);
                        } else {
                            edgecigar.push_back(cur_cigar);
                        }
                        DEBUG("CIGAR: " << a_i << " " << n << " " << edgeblocks[cur_block] << " " << edgecigar[edgecigar.size() - 1] << " " << i << " " << cigar.size());
                        a_i = edgeblocks[cur_block];
                        cur_cigar = "";
                        cur_block ++;
                        if (cur_block > edgeblocks.size()) {
                            WARN("CIGAR: Blocks ended! Something wrong with CIGAR alignment");
                            break;
                        } 
                    }
                    a_i += n;
                }
                if (c != 'I'){  
                    r_i += n;   
                    cur_end_pos += n;
                }
                cur_cigar += std::to_string(n) + c;
                cur_num = "";
            }
        }
        if (cur_cigar != "") {
            edgecigar.push_back(cur_cigar);
            edge_initial_ranges.push_back(Range(cur_start_pos, cur_end_pos));
            DEBUG("CIGAR: bounds  " << cur_start_pos << " " << cur_end_pos << " " << start << " " << r_i);
        }
    }

    void getMappedString(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappingpath, string &aligned, std::vector<size_t> &edgeblocks) {
        for (size_t i = 0; i < mappingpath.size(); ++ i) {
            EdgeId edgeid = mappingpath.edge_at(i);
            omnigraph::MappingRange mapping = mappingpath.mapping_at(i);
            size_t mapping_start = mapping.mapped_range.start_pos;
            size_t mapping_end = mapping.mapped_range.end_pos + gp_.g.k();
            if (i > 0){
                mapping_start = 0;
            }
            if (i < mappingpath.size() - 1) {
                mapping_end = gp_.g.length(edgeid);
            }
            string tmp = gp_.g.EdgeNucls(edgeid).str();
            string to_add = tmp.substr(mapping_start, mapping_end - mapping_start);
            aligned += to_add;
            edgeblocks.push_back(aligned.size());
        }
        return;
    }

    void getMappingOnRead(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappingpath, size_t &start, size_t &end) {
        start = mappingpath.mapping_at(0).initial_range.start_pos;
        end = mappingpath.mapping_at(mappingpath.size() - 1).initial_range.end_pos + gp_.g.k();
        return;
    }

    std::string getSubRead(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappingpath, const io::SingleRead &read) {
        size_t start;
        size_t end;
        getMappingOnRead(mappingpath, start, end);
        std::string readStr = read.sequence().str();
        return readStr.substr(start, end - start);
    }

    void ToGPA(const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &aligned_mappings, const io::SingleRead &read) {
        int nameIndex = 0;
        std::string res = "";
        for (const auto &mappingpath : aligned_mappings){
            string prev = "";
            string subread = getSubRead(mappingpath, read);
            string alignment;
            std::vector<size_t> edgeblocks;
            getMappedString(mappingpath, alignment, edgeblocks);
            std::vector<string>  edgecigar;
            std::vector<Range> edge_initial_ranges;
            int score;
            size_t start;
            size_t end;
            getMappingOnRead(mappingpath, start, end);
            getByEdgeCIGAR(subread, alignment, edgeblocks, start, edgecigar, edge_initial_ranges, score);

            for (size_t i = 0; i < mappingpath.size(); ++ i) {
                EdgeId edgeid = mappingpath.edge_at(i);
                omnigraph::MappingRange mapping = mappingpath.mapping_at(i);
                size_t mapping_start = mapping.mapped_range.start_pos;
                size_t mapping_end = mapping.mapped_range.end_pos + gp_.g.k();
                if (i > 0){
                    mapping_start = 0;
                }
                if (i < mappingpath.size() - 1) {
                    mapping_end = gp_.g.length(edgeid);
                }
                map<string, string> line = {{"Ind", "A"}, {"Name", ""}, {"ReadName", read.name()}, {"StartR", ""}, {"LenR", ""}, {"DirR", ""}
                                                                      , {"EdgeId", ""}, {"StartE", ""}, {"LenE", ""}, {"DirE", ""}
                                                                      , {"CIGAR", ""}, {"Prev", ""} , {"Next", ""}};
                
                nameIndex ++;                
                line["Name"] = read.name() + "_" + std::to_string(nameIndex);
                
                line["StartR"] = std::to_string(edge_initial_ranges[i].start_pos); 
                line["LenR"] = std::to_string(edge_initial_ranges[i].end_pos - 1 - edge_initial_ranges[i].start_pos); 
                line["DirR"] = "?"; // TODO


                line["EdgeId"] = std::to_string(edgeid.int_id());
                line["StartE"] = std::to_string(mapping_start);
                line["LenE"] = std::to_string(mapping_end - mapping_start);
                line["DirE"] = "?"; // TODO

                line["CIGAR"] = "";//edgecigar[i];

                if (i > 0){
                    line["Prev"] = prev;
                } else {
                    line["Prev"] = "-";
                }
                prev = line["Name"]; 
                if (i < mappingpath.size() - 1){
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
    
    bool IsCanonical(EdgeId e) const {
        return e <= gp_.g.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : gp_.g.conjugate(e);
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
            //ToGPA(aligned_mappings, read);
            string pathStr = "";
            for (const auto &mappingpath : aligned_mappings){
                for (const auto &edgeid: mappingpath.simple_path()) {
                    VertexId v1 = gp_.g.EdgeStart(edgeid);
                    VertexId v2 = gp_.g.EdgeEnd(edgeid);
                    pathStr += std::to_string(Canonical(edgeid).int_id()) + " (" + std::to_string(v1.int_id()) + "," + std::to_string(v2.int_id()) + ") ";
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
                size_t mapping_start = 0;
                size_t mapping_end = 0;
                string cur_path = "";
                string cur_path_len = "";
                string cur_substr = "";
                string str = "";
                for (size_t i = 0; i < path.size(); ++ i) {
                    EdgeId edgeid = path.edge_at(i);
                    omnigraph::MappingRange mapping = path.mapping_at(i);
                    mapping_start = mapping.mapped_range.start_pos;
                    mapping_end = mapping.mapped_range.end_pos + gp_.g.k();
                    if (i > 0){
                        mapping_start = 0;
                    }
                    if (i < path.size() - 1) {
                        mapping_end = gp_.g.length(edgeid);
                    }
                    cur_path += std::to_string(Canonical(edgeid).int_id()) + " (" + std::to_string(mapping_start) + "," + std::to_string(mapping_end) + ") ["
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
                pathStr += cur_path + "; ";
                subStr += cur_substr + "\n";
                int d = max((int) read.sequence().size(), 20);
                edlib::EdlibAlignResult result = edlib::edlibAlign(read.sequence().str().c_str(), (int) read.sequence().size(), str.c_str(), (int) str.size()
                                                   , edlib::edlibNewAlignConfig(d, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_DISTANCE,
                                                                         NULL, 0));
                int score = pacbio::STRING_DIST_INF;
                if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
                    score = result.editDistance;
                }
                edlib::edlibFreeAlignResult(result);
                if (seq_end - seq_start > max_len){
                    max_len = seq_end - seq_start;
                    max_str = read.name() + "\t" + std::to_string(seq_start) + "\t" + std::to_string(seq_end) + "\t"  
                                                    + std::to_string(read.sequence().size())+  "\t" + cur_path + "\t" + cur_path_len + "\t"+ std::to_string(score) + "\t" + str + "\n";
                }
                sum_str += read.name() + "\t" + std::to_string(seq_start) + "\t" + std::to_string(seq_end) + "\t" 
                                                 + std::to_string(read.sequence().size())+  "\t" + cur_path + "\t" + cur_path_len + "\t" + std::to_string(score) + "\t" + str + "\n";
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

void Launch(size_t K, const string &saves_path, const string &sequence_fasta, const string &mapper_type, const string &output_file, int threads) {
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

#pragma omp parallel num_threads(threads)
#pragma omp for
        for (size_t i =0 ; i < wrappedreads.size(); ++i) {
            aligner.AlignRead(wrappedreads[i]);
        }
}
}

int main(int argc, char **argv) {
    if (argc < 6) {
        cout << "Usage: longreads_aligner <K>"
             << " <saves path> <sequences file (fasta/fastq)> <mapper type: {pacbio, nanopore, default} > <ouput-prefix> {threads-num(16)}" << endl;
        exit(1);
    }
    int threads = 16;
    if (argc == 7) {
        threads = std::stoi(argv[6]);
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
    string output_file = argv[5];
    debruijn_graph::Launch(K, saves_path, sequence_file, mapper_type, output_file, threads);
    return 0;
}
