//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "mapping_printer.hpp"

#include "edlib/edlib.h"

namespace sensitive_aligner {

using namespace std;

string getStrId(const EdgeId &e, const debruijn_graph::ConjugateDeBruijnGraph &g_) {
    if (e.int_id() < g_.conjugate(e).int_id()) {
        return to_string(e.int_id()) + "+";
    } else {
        return to_string(e.int_id()) + "-";
    }
}

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
            path_str += getStrId(mappingpath[i], g_) + ",";
            path_len_str += to_string(mapping_end - mapping_start) + ",";
            path_seq_str += g_.EdgeNucls(mappingpath[i]).Subseq(mapping_start, mapping_end).str();
        }
        seq_starts += to_string(aligned_mappings.read_ranges[j].path_start.seq_pos) + ",";
        seq_ends += to_string(aligned_mappings.read_ranges[j].path_end.seq_pos) + ",";
        edge_starts += to_string(aligned_mappings.read_ranges[j].path_start.edge_pos) + ",";
        edge_ends += to_string(aligned_mappings.read_ranges[j].path_end.edge_pos) + ",";
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
            bwa_path_str += getStrId(edgeid, g_) + " (" + to_string(mapping.mapped_range.start_pos) + ","
                            + to_string(mapping.mapped_range.end_pos) + ") ["
                            + to_string(mapping.initial_range.start_pos) + ","
                            + to_string(mapping.initial_range.end_pos) + "], ";
        }
        bwa_path_str += ";";
    }
    string str = read.name() + "\t" + seq_starts + "\t"
                 + seq_ends + "\t"
                 + edge_starts + "\t"
                 + edge_ends + "\t"
                 + to_string(read.sequence().size()) +  "\t"
                 + path_str + "\t" + path_len_str + "\t" + bwa_path_str + "\t" + path_seq_str + "\n";
    DEBUG("Read " << read.name() << " aligned and length=" << read.sequence().size());
    DEBUG("Read " << read.name() << ". Paths with ends: " << path_str );
    #pragma omp critical
    {
        output_file_ << str;
    }
}



string MappingPrinterGPA::Print(map<string, string> &line) {
    vector<string> v = {"Ind", "Name", "ReadName", "StartR", "LenR", "DirR", "EdgeId", "StartE", "LenE", "DirE", "CIGAR", "Prev", "Next"};
    string outStr = "";
    for (const auto &it : v) {
        outStr += line[it] + "\t";
    }
    return outStr;
}


string MappingPrinterGPA::getCigar(const string &read, const string &aligned) {
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

void MappingPrinterGPA::getEdgeCigar(const string &subread, const string &path_seq, const vector<size_t> &edgeblocks,
                                     vector<string> &edgecigar, vector<Range> &edgeranges) {
    edgecigar.clear();
    edgeranges.clear();
    string cigar = getCigar(subread, path_seq);
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


void MappingPrinterGPA::getPath(const vector<debruijn_graph::EdgeId> &path,
                                const PathRange &path_range,
                                string &aligned, vector<size_t> &edgeblocks) {
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


string MappingPrinterGPA::getSubread(const Sequence &read, const PathRange &path_range) {
    return read.Subseq(path_range.path_start.seq_pos, path_range.path_end.seq_pos).str();
}

string MappingPrinterGPA::formGPAOutput(const io::SingleRead &read,
                                        const vector<debruijn_graph::EdgeId> &path,
                                        const vector<string> &edgecigar,
                                        const vector<Range> &edgeranges,
                                        int &nameIndex, const PathRange &path_range) {
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
            {"EdgeId", getStrId(edgeid, g_)},
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
    for (size_t i = 0; i < aligned_mappings.main_storage.size(); ++ i) {
        auto &path = aligned_mappings.main_storage[i];
        auto &path_range = aligned_mappings.read_ranges[i];

        string path_seq;
        vector<size_t> path_edgeblocks;
        getPath(path, path_range, path_seq, path_edgeblocks);

        string subread = getSubread(read.sequence(), path_range);

        vector<string> path_edgecigar;
        vector<Range> path_edgeranges;
        getEdgeCigar(subread, path_seq, path_edgeblocks, path_edgecigar, path_edgeranges);

        string path_line = formGPAOutput(read, path, path_edgecigar, path_edgeranges, nameIndex, path_range);
        #pragma omp critical
        {
            output_file_ << path_line;
        }
    }
}


} // namespace sensitive_aligner