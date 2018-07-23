#include <fstream>

#include "edlib/edlib.h"

#include "mapping_printer.hpp"

namespace sensitive_aligner {

void MappingPrinterTSV::SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) {
    string path_str = "";
    string path_len_str = "";
    string path_seq_str = "";
    string seq_starts = "";
    string seq_ends = "";
    string edge_starts = "";
    string edge_ends = "";
    for (size_t j = 0; j < aligned_mappings.main_storage.size(); ++ j) {
        auto &mappingpath = aligned_mappings.main_storage[j];
        for (size_t i = 0; i < mappingpath.size(); ++ i) {
            size_t mapping_start = i == 0 ? aligned_mappings.read_ranges[j].path_start.edge_pos : 0;
            size_t mapping_end = i == mappingpath.size() - 1 ?  aligned_mappings.read_ranges[j].path_end.edge_pos : g_.length(mappingpath[i]);
            path_str += std::to_string(mappingpath[i].int_id()) + ",";
            path_len_str += std::to_string(mapping_end - mapping_start) + ",";
            path_seq_str += g_.EdgeNucls(mappingpath[i]).Subseq(mapping_start, mapping_end).str();
        }
        seq_starts += std::to_string(aligned_mappings.read_ranges[j].path_start.seq_pos) + ",";
        seq_ends += std::to_string(aligned_mappings.read_ranges[j].path_end.seq_pos) + ",";
        edge_starts += std::to_string(aligned_mappings.read_ranges[j].path_start.edge_pos) + ",";
        edge_ends += std::to_string(aligned_mappings.read_ranges[j].path_end.edge_pos) + ",";
        path_str += ";";
        path_len_str += ";";
        path_seq_str += ";";
    }
    DEBUG("Paths: " << path_str);
    string bwa_path_str = "";
    for (const auto &path : aligned_mappings.bwa_paths) {
        for (size_t i = 0; i < path.size(); ++ i) {
            EdgeId edgeid = path.edge_at(i);
            omnigraph::MappingRange mapping = path.mapping_at(i);
            bwa_path_str += std::to_string(edgeid.int_id()) + " (" + std::to_string(mapping.mapped_range.start_pos) + ","
                            + std::to_string(mapping.mapped_range.end_pos) + ") ["
                            + std::to_string(mapping.initial_range.start_pos) + ","
                            + std::to_string(mapping.initial_range.end_pos) + "], ";
        }
        bwa_path_str += ";";
    }
    string str = read.name() + "\t" + seq_starts + "\t"
                 + seq_ends + "\t"
                 + edge_starts + "\t"
                 + edge_ends + "\t"
                 + std::to_string(read.sequence().size()) +  "\t" 
                 + path_str + "\t" + path_len_str + "\t" + bwa_path_str + "\t" + path_seq_str + "\n";
    DEBUG("Read " << read.name() << " aligned and length=" << read.sequence().size());
    DEBUG("Read " << read.name() << ". Paths with ends: " << path_str );
    #pragma omp critical
    {
        output_file_ << str;
    }
}



std::string MappingPrinterGPA::Print(map<string, string> &line) {
    std::vector<string> v = {"Ind", "Name", "ReadName", "StartR", "LenR", "DirR", "EdgeId", "StartE", "LenE", "DirE", "CIGAR", "Prev", "Next"};
    string outStr = "";
    for (const auto &it : v) {
        outStr += line[it] + "\t";
    }
    return outStr;
}


void MappingPrinterGPA::CIGAR(std::string &read, std::string aligned, std::string &cigar, int &score) {
    int d = max((int) read.size(), 20);
    edlib::EdlibAlignResult result = edlib::edlibAlign(aligned.c_str(), (int) aligned.size(), read.c_str(), (int) read.size()
                                     , edlib::edlibNewAlignConfig(d, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_PATH,
                                             NULL, 0));
    cigar = "";
    score = std::numeric_limits<int>::max();
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
            if (c == '=' || c == 'I' || c == 'X' || c == 'M') {
                len_a += n;
            }
            if (c != 'I') {
                len_r += n;
            }
            cur_num = "";
            n = 0;
        }
    }
    DEBUG("CIGAR: " << len_a << " " << aligned.size()  << " " << len_r << " " << read.size());
}

void MappingPrinterGPA::DivideByEdgeCIGAR(string &read, string &aligned, std::vector<size_t> &edgeblocks, size_t start,
                                       std::vector<string> &edgecigar, std::vector<Range> &edge_initial_ranges, int &score) {
    std::string cigar;
    CIGAR(read, aligned, cigar, score);
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
            if (c == '=' || c == 'I' || c == 'X' || c == 'M') {
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
            if (c != 'I') {
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

void MappingPrinterGPA::MappedString(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappingpath, string &aligned, std::vector<size_t> &edgeblocks) {
    for (size_t i = 0; i < mappingpath.size(); ++ i) {
        EdgeId edgeid = mappingpath.edge_at(i);
        omnigraph::MappingRange mapping = mappingpath.mapping_at(i);
        size_t mapping_start = mapping.mapped_range.start_pos;
        size_t mapping_end = mapping.mapped_range.end_pos + g_.k();
        if (i > 0) {
            mapping_start = 0;
        }
        if (i < mappingpath.size() - 1) {
            mapping_end = g_.length(edgeid);
        }
        string tmp = g_.EdgeNucls(edgeid).str();
        string to_add = tmp.substr(mapping_start, mapping_end - mapping_start);
        aligned += to_add;
        edgeblocks.push_back(aligned.size());
    }
    return;
}

void MappingPrinterGPA::MappingOnRead(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappingpath, size_t &start, size_t &end) {
    start = mappingpath.mapping_at(0).initial_range.start_pos;
    end = mappingpath.mapping_at(mappingpath.size() - 1).initial_range.end_pos + g_.k();
    return;
}

std::string MappingPrinterGPA::SubRead(const omnigraph::MappingPath<debruijn_graph::EdgeId> &mappingpath, const io::SingleRead &read) {
    size_t start;
    size_t end;
    MappingOnRead(mappingpath, start, end);
    std::string readStr = read.sequence().str();
    return readStr.substr(start, end - start);
}

void MappingPrinterGPA::SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) {
    int nameIndex = 0;
    std::string res = "";
    for (const auto &mappingpath : aligned_mappings.bwa_paths) {
        string prev = "";
        string subread = SubRead(mappingpath, read);
        string alignment;
        std::vector<size_t> edgeblocks;
        MappedString(mappingpath, alignment, edgeblocks);
        std::vector<string>  edgecigar;
        std::vector<Range> edge_initial_ranges;
        int score;
        size_t start;
        size_t end;
        MappingOnRead(mappingpath, start, end);
        DivideByEdgeCIGAR(subread, alignment, edgeblocks, start, edgecigar, edge_initial_ranges, score);

        for (size_t i = 0; i < mappingpath.size(); ++ i) {
            EdgeId edgeid = mappingpath.edge_at(i);
            omnigraph::MappingRange mapping = mappingpath.mapping_at(i);
            size_t mapping_start = mapping.mapped_range.start_pos;
            size_t mapping_end = mapping.mapped_range.end_pos + g_.k();
            if (i > 0) {
                mapping_start = 0;
            }
            if (i < mappingpath.size() - 1) {
                mapping_end = g_.length(edgeid);
            }
            map<string, string> line = {{"Ind", "A"}, {"Name", ""}, {"ReadName", read.name()}, {"StartR", ""}, {"LenR", ""}, {"DirR", ""}
                , {"EdgeId", ""}, {"StartE", ""}, {"LenE", ""}, {"DirE", ""}
                , {"CIGAR", ""}, {"Prev", ""} , {"Next", ""}
            };

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

            if (i > 0) {
                line["Prev"] = prev;
            } else {
                line["Prev"] = "-";
            }
            prev = line["Name"];
            if (i < mappingpath.size() - 1) {
                line["Next"] = read.name() + "_" + std::to_string(nameIndex + 1);
            } else {
                line["Next"] = "-";
            }
            res += Print(line) + "\n";
        }
    }
    #pragma omp critical
    {
        output_file_ << res;
    }
}

} // namespace sensitive_aligner