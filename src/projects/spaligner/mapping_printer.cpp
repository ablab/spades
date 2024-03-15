//***************************************************************************
//* Copyright (c) 2018-2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "mapping_printer.hpp"

#include "edlib/edlib.h"

#include <sstream>

namespace sensitive_aligner {

using namespace std;

string MappingPrinter::StrId(const EdgeId &e) const {
    string id_str = edge_namer_.EdgeOrientationString(e);
    return id_str;
}

void MappingPrinterTSV::SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) {
    stringstream path_ss;
    stringstream path_len_ss;
    stringstream path_seq_ss;
    stringstream seq_starts;
    stringstream seq_ends;
    stringstream edge_starts;
    stringstream edge_ends;
    for (size_t j = 0; j < aligned_mappings.edge_paths.size(); ++ j) {
        auto &mappingpath = aligned_mappings.edge_paths[j];
        for (size_t i = 0; i < mappingpath.size(); ++ i) {
            size_t mapping_start = i == 0 ? aligned_mappings.read_ranges[j].path_start.edge_pos : 0;
            size_t mapping_end = i == mappingpath.size() - 1 ?  aligned_mappings.read_ranges[j].path_end.edge_pos : g_.length(mappingpath[i]);
            string delim = i == mappingpath.size() - 1 ? "" : ",";
            path_ss << StrId(mappingpath[i]) << delim;
            path_len_ss << mapping_end - mapping_start << delim;
            path_seq_ss << g_.EdgeNucls(mappingpath[i]).Subseq(mapping_start, mapping_end).str();
        }
        string delim = j == aligned_mappings.edge_paths.size() - 1 ? "" : ",";
        seq_starts << aligned_mappings.read_ranges[j].path_start.seq_pos << delim;
        seq_ends << aligned_mappings.read_ranges[j].path_end.seq_pos << delim;
        edge_starts << aligned_mappings.read_ranges[j].path_start.edge_pos << delim;
        edge_ends << aligned_mappings.read_ranges[j].path_end.edge_pos << delim;
        string s_delim = j == aligned_mappings.edge_paths.size() - 1 ? "" : ";";
        path_ss << s_delim;
        path_len_ss << s_delim;
        path_seq_ss << s_delim;
    }
    string str = read.name() + "\t" + seq_starts.str() + "\t"
                 + seq_ends.str() + "\t"
                 + edge_starts.str() + "\t"
                 + edge_ends.str() + "\t"
                 + to_string(read.sequence().size()) +  "\t"
                 + path_ss.str() + "\t" + path_len_ss.str() + "\t" + path_seq_ss.str() + "\n";
    DEBUG("Read " << read.name() << " aligned and length=" << read.sequence().size());
    #pragma omp critical
    {
        output_file_ << str;
    }
}

void MappingPrinterFasta::SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) {
    string str = "";
    for (size_t j = 0; j < aligned_mappings.edge_paths.size(); ++ j) {
        auto &mappingpath = aligned_mappings.edge_paths[j];
        string path_str = "";
        string path_seq_str = "";
        for (size_t i = 0; i < mappingpath.size(); ++ i) {
            size_t mapping_start = i == 0 ? aligned_mappings.read_ranges[j].path_start.edge_pos : 0;
            size_t mapping_end = i == mappingpath.size() - 1 ?  aligned_mappings.read_ranges[j].path_end.edge_pos : g_.length(mappingpath[i]);
            string delim = i == mappingpath.size() - 1 ? "" : "_";
            path_str += StrId(mappingpath[i]) + delim;
            path_seq_str += g_.EdgeNucls(mappingpath[i]).Subseq(mapping_start, mapping_end).str();
        }
        str += ">" + read.name() + "|Edges=" + path_str
                                 + "|start_g=" + to_string(aligned_mappings.read_ranges[j].path_start.edge_pos)
                                 + "|end_g=" + to_string(aligned_mappings.read_ranges[j].path_end.edge_pos)
                                 + "|start_s=" + to_string(aligned_mappings.read_ranges[j].path_start.seq_pos)
                                 + "|end_s=" + to_string(aligned_mappings.read_ranges[j].path_end.seq_pos)
                                 + "\n" + path_seq_str + "\n";
    }
    #pragma omp critical
    {
        output_file_ << str;
    }
}

string MappingPrinterGPA::Print(map<string, string> &line) const {
    vector<string> v = {"Ind", "Name", "ReadName", "StartR", "LenR", "DirR", "EdgeId", "StartE", "LenE", "DirE", "CIGAR", "Prev", "Next"};
    string outStr = "";
    for (const auto &it : v) {
        outStr += line[it] + "\t";
    }
    return outStr;
}


string MappingPrinterGPA::FormCigar(const string &read, const string &aligned) const {
    int d = max((int) read.size(), 20);
    edlib::EdlibAlignResult result = edlib::edlibAlign(aligned.c_str(), (int) aligned.size(), read.c_str(), (int) read.size()
                                     , edlib::edlibNewAlignConfig(d, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_PATH,
                                             NULL, 0));
    string cigar = "";
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        cigar = edlib::edlibAlignmentToCigar(result.alignment, result.alignmentLength, edlib::EDLIB_CIGAR_EXTENDED);
    }
    edlib::edlibFreeAlignResult(result);
    return cigar;
}

void MappingPrinterGPA::FormEdgeCigar(const string &subread, const string &path_seq,
                                      const vector<size_t> &edgeblocks,
                                      vector<string> &edgecigar,
                                      vector<Range> &edgeranges) const {
    edgecigar.clear();
    edgeranges.clear();
    string cigar = FormCigar(subread, path_seq);
    string cur_num = "";
    int cur_block = 0;
    int seq_start = 0;
    int seq_end = 0;
    string cur_cigar = "";
    int a_i = 0;
    for (size_t i = 0; i < cigar.size(); ++ i) {
        if (isdigit(cigar[i])) {
            cur_num += cigar[i];
        } else {
            int n = stoi(cur_num);
            cur_num = "";
            char c = cigar[i];

            if (c == '=' || c == 'I' || c == 'X' || c == 'M') {
                while (a_i + n > (int) edgeblocks[cur_block]) {
                    DEBUG("CIGAR: " << n << c);
                    n -= (int) (edgeblocks[cur_block] - a_i);
                    if (c != 'I') {
                        seq_end += (int) (edgeblocks[cur_block] - a_i);
                    }
                    edgeranges.push_back(Range(seq_start, seq_end));
                    edgecigar.push_back(cur_cigar + to_string((int)edgeblocks[cur_block] - a_i) + c);
                    DEBUG("CIGAR: " << a_i << " " << n << " " << edgeblocks[cur_block] << " " << edgecigar[edgecigar.size() - 1] << " " << i << " " << cigar.size());
                    a_i = (int) edgeblocks[cur_block];
                    seq_start = seq_end;
                    cur_cigar = "";
                    cur_block ++;
                    if (cur_block > (int) edgeblocks.size()) {
                        WARN("CIGAR: Blocks ended! Something wrong with CIGAR alignment");
                        break;
                    }
                }
                a_i += n;
            }
            if (c != 'I') {
                seq_end += n;
            }
            cur_cigar += to_string(n) + c;
        }
    }
    if (cur_cigar != "") {
        edgecigar.push_back(cur_cigar);
        edgeranges.push_back(Range(seq_start, seq_end));
        DEBUG("CIGAR: bounds  " << seq_start << " " << seq_end);
    }

}


void MappingPrinterGPA::GeneratePath(const vector<debruijn_graph::EdgeId> &path,
                                const PathRange &path_range,
                                string &aligned, vector<size_t> &edgeblocks) const {
    aligned = "";
    edgeblocks.clear();
    for (size_t i = 0; i < path.size(); ++ i) {
        EdgeId edgeid = path[i];
        size_t mapping_start = 0;
        size_t mapping_end = g_.length(edgeid);
        if (i == 0) {
            mapping_start = path_range.path_start.edge_pos;
        }
        if (i == path.size() - 1) {
            mapping_end = path_range.path_end.edge_pos;
        }
        aligned += g_.EdgeNucls(edgeid).Subseq(mapping_start, mapping_end).str();
        edgeblocks.push_back(aligned.size());
    }
    return;
}


string MappingPrinterGPA::GenerateSubread(const Sequence &read, const PathRange &path_range) const {
    return read.Subseq(path_range.path_start.seq_pos, path_range.path_end.seq_pos).str();
}

string MappingPrinterGPA::FormGPAOutput(const io::SingleRead &read,
                                        const vector<debruijn_graph::EdgeId> &path,
                                        const vector<string> &edgecigar,
                                        const vector<Range> &edgeranges,
                                        int &nameIndex, const PathRange &path_range) const {
    string res;
    string prev = "-";
    string next = "-";
    for (size_t i = 0; i < path.size(); ++ i) {
        EdgeId edgeid = path[i];
        size_t start = i == 0 ? path_range.path_start.edge_pos : 0;
        size_t end = i + 1 == path.size() ? path_range.path_end.edge_pos : g_.length(edgeid);
        next = i + 1 == path.size() ? "-" : read.name() + "_" + to_string(nameIndex + 1);
        map<string, string> line = {{"Ind", "A"},
            {"Name", read.name() + "_" + to_string(nameIndex)},
            {"ReadName", read.name()},
            {"StartR", to_string(path_range.path_start.seq_pos + edgeranges[i].start_pos)},
            {"LenR", to_string(edgeranges[i].end_pos - edgeranges[i].start_pos)},
            {"DirR", "+"},
            {"EdgeId", StrId(edgeid)},
            {"StartE", to_string(start)},
            {"LenE", to_string(end - start)},
            {"DirE", "+"},
            {"CIGAR", edgecigar[i]},
            {"Prev", prev},
            {"Next", next}
        };
        prev = read.name() + "_" + to_string(nameIndex);
        ++ nameIndex;
        res += Print(line) + "\n";
    }

    return res;

}

void MappingPrinterGPA::SaveMapping(const sensitive_aligner::OneReadMapping &aligned_mappings, const io::SingleRead &read) {
    int nameIndex = 0;
    for (size_t i = 0; i < aligned_mappings.edge_paths.size(); ++ i) {
        auto &path = aligned_mappings.edge_paths[i];
        auto &path_range = aligned_mappings.read_ranges[i];

        string path_seq;
        vector<size_t> path_edgeblocks;
        GeneratePath(path, path_range, path_seq, path_edgeblocks);

        string subread = GenerateSubread(read.sequence(), path_range);

        vector<string> path_edgecigar;
        vector<Range> path_edgeranges;
        FormEdgeCigar(subread, path_seq, path_edgeblocks, path_edgecigar, path_edgeranges);

        string path_line = FormGPAOutput(read, path, path_edgecigar, path_edgeranges, nameIndex, path_range);
        #pragma omp critical
        {
            output_file_ << path_line;
        }
    }
}


} // namespace sensitive_aligner