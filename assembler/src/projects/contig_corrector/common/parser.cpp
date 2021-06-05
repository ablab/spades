//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "common.hpp"

using namespace std;

namespace {

bool is_next_contig(string const & s) {
    return s[0] == '>';
}

pair<string, string> getNameAndInfo(string const & s) {
    // s looks like:
    // >tig00000199 len=181991 reads=225 class=contig suggestRepeat=no suggestBubble=no suggestCircular=no
    pair<string, string> name_and_info;
    auto name_end = s.find(' ', 1);
    if (name_end == string::npos)
        name_end = s.size();
    name_and_info.first = s.substr(1, name_end - 1); // skip '>';

    if (name_end + 1 < s.size())
        name_and_info.second = s.substr(name_end + 1);
    return name_and_info;
}

} // namespace

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
            std::tie(current_seq.name, current_seq.info) = getNameAndInfo(current_line);
            continue;
        }

        current_seq.seq.append(current_line);
    }

    if (!current_seq.name.empty())
        contigs.push_back(move(current_seq));
    return contigs;
}
