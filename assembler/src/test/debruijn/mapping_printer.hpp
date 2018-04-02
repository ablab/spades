#include "edlib/edlib.h"

namespace debruijn_graph {

class MappingPrinter {
protected:
    const conj_graph_pack &gp_;
    std::string output_file_;

public:

    MappingPrinter(const conj_graph_pack &gp, const std::string &output_file)
                :gp_(gp), output_file_(output_file)
    {}

    virtual void SaveMapping(const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &aligned_mappings, const io::SingleRead &read) = 0;

    virtual ~MappingPrinter (){};

};

class MappingPrinterTSV: public MappingPrinter {

public:

    MappingPrinterTSV(const conj_graph_pack &gp, const std::string &output_file)
        :MappingPrinter(gp, output_file)
    {
        ofstream tsv_file;
        tsv_file.open(output_file_ + ".tsv", std::ofstream::out);
        //tsv_file << "read\tstart_pos\tend_pos\tread_length\tgraph_path\tedges_lengths\tmapping\ted\n";
        tsv_file.close();
    }

    virtual void SaveMapping(const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &aligned_mappings, const io::SingleRead &read) {
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
        string sum_str = "";
        string max_str = "";
        int max_len = 0;

        string cur_path_sum = "";
        string cur_path_len_sum = "";
        string cur_str_sum = "";
        int seq_start_sum = -1;
        int seq_end_sum = 0;

        for (const auto &path : aligned_mappings){
            size_t mapping_start = 0;
            size_t mapping_end = 0;
            string cur_str = "";
            string cur_path = "";
            string cur_path_len = "";
            int seq_start = -1;
            int seq_end = 0;
            for (size_t i = 0; i < path.size(); ++ i) {
                EdgeId edgeid = path.edge_at(i);
                omnigraph::MappingRange mapping = path.mapping_at(i);
                mapping_start = mapping.mapped_range.start_pos;
                mapping_end = mapping.mapped_range.end_pos + gp_.g.k();
                //initial_start = mapping.initial_range.start_pos;
                ///initial_end = mapping.initial_range.end_pos + gp_.g.k();
                if (i > 0){
                    mapping_start = 0;
                }
                if (i < path.size() - 1) {
                    mapping_end = gp_.g.length(edgeid);
                    //initial_end = mapping.initial_range.end_pos;
                }
                cur_path += std::to_string(edgeid.int_id()) + " (" + std::to_string(mapping_start) + "," + std::to_string(mapping_end) + ") ["
                           + std::to_string(mapping.initial_range.start_pos) + "," + std::to_string(mapping.initial_range.end_pos) + "], ";
                cur_path_len += std::to_string(mapping_end - mapping_start) + ",";
                if (seq_start < 0){
                    seq_start = (int) mapping.initial_range.start_pos;
                }
                seq_end = (int) mapping.initial_range.end_pos;

                cur_path_sum += std::to_string(edgeid.int_id()) + " (" + std::to_string(mapping_start) + "," + std::to_string(mapping_end) + ") ["
                           + std::to_string(mapping.initial_range.start_pos) + "," + std::to_string(mapping.initial_range.end_pos) + "], ";
                cur_path_len_sum += std::to_string(mapping_end - mapping_start) + ",";
                if (seq_start_sum < 0){
                    seq_start_sum = (int) mapping.initial_range.start_pos;
                }
                seq_end_sum = (int) mapping.initial_range.end_pos;

                string tmp = gp_.g.EdgeNucls(edgeid).str();
                string to_add = tmp.substr(mapping_start, mapping_end - mapping_start);
                cur_str += to_add;
            }
            cur_str_sum += cur_str;
            if (seq_end - seq_start > max_len) {
                max_len = seq_end - seq_start;
                edlib::EdlibAlignResult result = edlib::edlibAlign(cur_str.c_str(), (int) cur_str.size(), read.sequence().str().c_str(), (int) read.sequence().size()
                                                   , edlib::edlibNewAlignConfig(-1, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_DISTANCE,
                                                                         NULL, 0));
                //if (result.editDistance >= 0 && max_len >= 1200) {
                    max_str = read.name() + "\t" + std::to_string(seq_start) + "\t" + std::to_string(seq_end) + "\t" 
                                                     + std::to_string(read.sequence().size())+  "\t" + cur_path + "\t" + cur_path_len + "\t" + cur_str + "\t" 
                                                     + std::to_string(result.editDistance) + "\n";
                //}
                edlib::edlibFreeAlignResult(result);
            }
        }
        sum_str = read.name() + "\t" + std::to_string(seq_start_sum) + "\t" + std::to_string(seq_end_sum) + "\t" 
                                             + std::to_string(read.sequence().size())+  "\t" + cur_path_sum + "\t" + cur_path_len_sum + "\t" + cur_str_sum + "\n";
        INFO("Read " << read.name() << " aligned and length=" << read.sequence().size());
        INFO("Read " << read.name() << ". Paths with ends: " << cur_path_sum );
        //INFO("Seq subs: " << subStr);
#pragma omp critical
        {
            //if (max_len >= 1200) {
                ofstream myfile;
                myfile.open(output_file_ + ".tsv", std::ofstream::out | std::ofstream::app);
                myfile << max_str;
                myfile.close();
            //}
        }
    }
    

};

class MappingPrinterGPA : public MappingPrinter {
public:
    MappingPrinterGPA(const conj_graph_pack &gp, const std::string &output_file)
        :MappingPrinter(gp, output_file)
    {
        ofstream gpa_file;
        gpa_file.open(output_file_ + ".gpa", std::ofstream::out);
        gpa_file << "H\n";
        gpa_file.close();
    }

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

    virtual void SaveMapping(const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &aligned_mappings, const io::SingleRead &read) {
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


};

class MappingPrinterHub {
private:
    std::string output_file_;
    vector<MappingPrinter*> mapping_printers_;

public:

    MappingPrinterHub(const conj_graph_pack &gp, const std::string &output_file, const std::string formats)
                :output_file_(output_file) {
        if (formats.find("tsv") != std::string::npos) {
            mapping_printers_.push_back(new MappingPrinterTSV(gp, output_file_));
        }
        if (formats.find("gpa") != std::string::npos) {
            mapping_printers_.push_back(new MappingPrinterGPA(gp, output_file_));
        }
    }

    void SaveMapping(const std::vector<omnigraph::MappingPath<debruijn_graph::EdgeId> > &aligned_mappings, const io::SingleRead &read) {
        for (auto printer: mapping_printers_) {
            printer->SaveMapping(aligned_mappings, read);
        }
    }

    ~MappingPrinterHub() {
        for (auto printer: mapping_printers_) {
            delete printer;
        }
        mapping_printers_.clear();
    }
};

}