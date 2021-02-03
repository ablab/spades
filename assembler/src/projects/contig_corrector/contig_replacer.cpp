//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contig_replacer.hpp"
#include "common.hpp"
#include "replacer.hpp"

#include <string>
#include <algorithm>
#include <unordered_map>

using namespace std;

namespace contig_replacer {

struct gcfg {
    std::string contigs_file;
    std::string fragments_file;
    std::string output_file;
} cfg;


clipp::group GetCLI() {
  using namespace clipp;

  auto cli = (
      cfg.contigs_file << value("contigs file"),
      cfg.fragments_file << value("fragments file"),
      cfg.output_file << value("output file")
  );

  return cli;
}

namespace {

void WriteWithWidth(std::ostream & out, std::string const & seq, size_t width = 50) {
    for (size_t pos = 0; pos < seq.size(); pos += width)
        out << seq.substr(pos, width) << '\n';
}

unordered_map<string, vector<size_t>> GroupByName(vector<SeqString> const & seqs) {
    unordered_map<string, vector<size_t>> res;
    for (size_t i = 0; i < seqs.size(); ++i)
        res[seqs[i].name].push_back(i);
    return res;
}

size_t ExtractNumber(std::string const & s, size_t from) {
    auto len = 0;
    while (from + len < s.size() && isdigit(s[from + len]))
        ++len;
    return stoull(s.substr(from, len));
}

using Bounds = pair<size_t, size_t>;

Bounds ExtractBounds(std::string const & info) {
    Bounds res;
    auto start_name = "start_s="s;
    auto end_name = "end_s="s;
    res.first = ExtractNumber(info, info.rfind(start_name) + start_name.size());
    res.second = ExtractNumber(info, info.rfind(end_name) + end_name.size());
    return res;
}

vector<Bounds> ExtractBounds(vector<SeqString> const & seqs) {
    vector<Bounds> res;
    res.reserve(seqs.size());
    for (auto const & s : seqs)
        res.push_back(ExtractBounds(s.info));
    return res;
}

void SortByBounds(unordered_map<string, vector<size_t>> & groups, vector<Bounds> const & bounds) {
    for (auto & group : groups)
        std::sort(group.second.begin(), group.second.end(),
                 [&bounds](size_t lhs, size_t rhs) { return bounds[lhs] < bounds[rhs]; });
}

list<ReplaceInfo> GenReplaceInfo(vector<SeqString> & fragments, vector<Bounds> const & bounds , vector<size_t> const & indexes) {
    list<ReplaceInfo> res;
    for (auto i : indexes)
        res.emplace_back(move(fragments[i].seq), bounds[i].first, bounds[i].second);
    return res;
}

} //namespace

int main(int argc, char * argv[]) {
    START_BANNER("SPAdes standalone contig replacer");

    auto contigs = ReadContigs(cfg.contigs_file);
    auto fragments = ReadContigs(cfg.fragments_file);
    auto fragment_groups = GroupByName(fragments);
    auto bounds = ExtractBounds(fragments);
    SortByBounds(fragment_groups, bounds);

    for (auto & contig : contigs) {
        auto it = fragment_groups.find(contig.name);
        if (it == fragment_groups.end())
            continue;
            INFO("Processing " << contig.name << " with " << it->second.size() << " fragments");
        contig.seq = Replace(contig.seq, GenReplaceInfo(fragments, bounds, it->second));
    }

    ofstream contigs_output(cfg.output_file);
    VERIFY(contigs_output.is_open());

    for (auto const & contig : contigs) {
        if (contig.seq.empty())
            continue;
        contigs_output << '>' << contig.name << " len=" << contig.seq.size() << '\n';
        WriteWithWidth(contigs_output, contig.seq);
    }

    INFO("SPAdes standalone contig replacer finished");

    return 0;
}

} // namespace contig_replacer
