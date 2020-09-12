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

pair<string, size_t> getNameAndLen(string const & s) {
    // s looks like:
    // >tig00000199 len=181991 reads=225 class=contig suggestRepeat=no suggestBubble=no suggestCircular=no
    auto name_end = s.find(' ', 1);
    auto name = s.substr(1, name_end - 1); // skip '>';

    std::string pattern = "len=";
    auto len_start = s.find(pattern, name_end + 1);
    if (len_start == std::string::npos)
        return {move(name), 0};
    size_t len_end = len_start;
    while (len_end < s.size() && isdigit(s[len_end]))
        ++len_end;
    if (len_start == len_end)
        return {move(name), 0};

    auto len = stoull(s.substr(len_start, len_end - len_start));
    return {move(name), len};
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
