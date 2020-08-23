//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "common.hpp"

using namespace debruijn_graph;
using namespace path_extend;
using namespace std;

namespace {

void AddEdgeId(const Graph& graph, const std::string& current_id, const std::unordered_map<size_t, EdgeId>& edge_map,
                                   std::vector<EdgeId>& current_alternative, bool rc = false)
{
    size_t edgeid = (size_t) std::stoi(current_id);
    auto e_iter = edge_map.find(edgeid);
    if (e_iter != edge_map.end()) {
        if (!rc)
            current_alternative.push_back(e_iter->second);
        else
            current_alternative.push_back(graph.conjugate(e_iter->second));

    } else {
        FATAL_ERROR("Edge id " << edgeid << " was not found");
    }
}

} // namespace

PathWithEdgePostionsContainer ParseInputPaths(const std::string& transcript_paths_file, const Graph& graph) {
    PathWithEdgePostionsContainer result;
    std::unordered_map<size_t, EdgeId> edge_map;
    for (auto iter = graph.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        edge_map.emplace(graph.int_id(*iter), *iter);
    }

    std::ifstream in_stream(transcript_paths_file);
    if (!in_stream.is_open())
        FATAL_ERROR("Failed to open " << transcript_paths_file);

    std::string line;
    while (std::getline(in_stream, line)) {
        size_t alternative_count = 1;
        //INFO(line)
        bool current_int_is_edgeid = true;
        std::string current_int;
        std::vector<EdgeId> current_alternative;
        PathWithEdgePostions current_path;

        for (size_t i = 0; i < line.length(); ++i) {
            if (std::isdigit(line[i]) || line[i] == '-') {
                current_int += line[i];
            } else if (line[i] == '\'') {
                if (current_int_is_edgeid) {
                    AddEdgeId(graph, current_int, edge_map, current_alternative, true);
                } else {
                    ERROR("Wrong format");
                }
                current_int_is_edgeid = true;
                current_int = "";
            } else if (line[i] == ',') {
                if (!current_int.empty()) {
                    if (current_int_is_edgeid) {
                        AddEdgeId(graph, current_int, edge_map, current_alternative);
                    } else {
                        ERROR("Wrong format");
                    }
                    current_int = "";
                    current_int_is_edgeid = true;
                }
            } else if (line[i] == ':') {
                if (!current_int.empty()) {
                    if (current_int_is_edgeid) {
                        AddEdgeId(graph, current_int, edge_map, current_alternative);
                    } else {
                        int pos = std::stoi(current_int);
                        current_path.positions.push_back(pos);
                    }
                    current_int = "";
                    current_int_is_edgeid = true;
                }
                alternative_count *= current_alternative.size();
                current_path.edge_set.emplace_back(std::move(current_alternative));
                current_alternative = std::vector<EdgeId>();
            } else if (line[i] == '_') {
                if (!current_int.empty()) {
                    if (current_int_is_edgeid) {
                        AddEdgeId(graph, current_int, edge_map, current_alternative);
                    } else {
                        ERROR("Wrong format");
                    }
                }
                current_int = "";
                current_int_is_edgeid = false;
            } else {
                WARN("Unknown charachter " << line[i]);
            }
        }
        if (!current_int.empty()) {
            if (current_int_is_edgeid) {
                AddEdgeId(graph, current_int, edge_map, current_alternative);
            } else {
                current_path.positions.push_back(std::stoi(current_int));
            }
        }

        alternative_count *= current_alternative.size();
        current_path.edge_set.emplace_back(std::move(current_alternative));
        VERIFY(current_path.positions.empty() || current_path.edge_set.size() == current_path.positions.size());

        if (alternative_count > 300) {
            DEBUG("The following path " << alternative_count <<  "  exceeds the limit for the number of alternatives");
            DEBUG(line)
            result.emplace_back();
        } else {
            result.emplace_back(std::move(current_path));
        }
    }
    return result;
}

bool is_next_contig(string const & s) {
    return s[0] == '>';
}

pair<string, size_t> getNameAndLen(string const & s) {
    // s looks like:
    // >tig00000199 len=181991 reads=225 class=contig suggestRepeat=no suggestBubble=no suggestCircular=no
    auto name_end = s.find(' ', 1);
    auto name = s.substr(1, name_end - 1); // skip '>';
    auto len_end = s.find(' ', name_end + 1);
    auto len_start = name_end + 5; // skip ' len='
    auto len = stoull(s.substr(len_start, len_end - len_start));
    return {move(name), len};
}

std::vector<SeqString> ReadContigs(std::string const & contigs_file) {
    ifstream inp(contigs_file);
    if (!inp.is_open())
        throw "Cannot open " + contigs_file;
    string current_line;
    vector<SeqString> contigs;
    SeqString current_seq;
    while (getline(inp, current_line)) {
        if (current_line.empty())
            continue;
        
        if (is_next_contig(current_line)) {
            if (!current_seq.name.empty())
                contigs.push_back(move(current_seq));
            auto name_and_len = getNameAndLen(current_line);
            current_seq.name = name_and_len.first;
            current_seq.seq.reserve(name_and_len.second);
            continue;
        }

        current_seq.seq.append(current_line);
    }

    if (!current_seq.name.empty())
        contigs.push_back(move(current_seq));
    return contigs;
}
