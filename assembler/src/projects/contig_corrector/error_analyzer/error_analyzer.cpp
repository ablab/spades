#include "error_analyzer.hpp"
#include "helpers/common.hpp"
#include "helpers/aligner_output_reader.hpp"

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <numeric>

enum class Columns : size_t {
    ref_name           = 0,
    contig_name        = 1,
    ref_start_pos      = 2,
    ref_nucls          = 3,
    contig_nucls       = 4,
    contig_start_pos   = 5,
    TOTAL_COLUMNS_SIZE = 6
};

template<Columns el>
struct type_getter<Columns, el> {
    using type = std::conditional_t<
                        el == Columns::ref_start_pos ||
                        el == Columns::contig_start_pos,
                        unsigned long long,
                        std::string>;
};

namespace error_analyzer {

struct gcfg {
    std::string contigs;
    std::string snps;
} cfg;


clipp::group GetCLI() {
  using namespace clipp;

  auto cli = (
      cfg.contigs << value("file with contigs"),
      cfg.snps << value("used_snps.gz")
  );

  return cli;
}

namespace {

template<Columns ... columns>
Records<Columns, columns ...> ReadUsedSnps(std::string const & file) {
    std::ifstream inp(file);
    if (!inp.is_open())
        throw "Cannot open '" + file + "' file";
    Records<Columns, columns ...> records;
    RecordPusher<Columns, columns ...> pusher(records, [](...){ return true; });
    std::string line;
    while (GetNextNonemptyLine(inp, line))
        pusher.Push(line);
    return records;
}

std::unordered_map<std::string, size_t> IndexByName(std::vector<SeqString> const & seqs) {
    std::unordered_map<std::string, size_t> res;
    for (size_t i = 0; i < seqs.size(); ++i) {
        auto success = res.emplace(seqs[i].name, i);
        if (!success.second)
            throw "The file with contigs contains a sequence with the non-unique name '" + seqs[i].name + "'";
    }
    return res;
}

enum class StatType {
    on_uncorrected = 0,
    on_corrected = 1,
    on_bound = 2
};

struct Stat : std::array<size_t, 3> {
    using Base = std::array<size_t, 3>;

    Stat() : Base({0,0,0}) {};

    size_t & operator[](StatType type) {
        return Base::operator[](static_cast<size_t>(type));
    }
    
    size_t operator[](StatType type) const {
        return Base::operator[](static_cast<size_t>(type));
    }
};

StatType GetStatType(std::string const &seq, size_t start_pos, size_t len) {
    bool ContainCorrected = false;
    bool ContainUncorrected = false;

    // extend checking range to left by 1 if possible
    if (start_pos > 0) {
        --start_pos;
        ++len;
    }

    // extend checking range to right by 1 if possible
    if (start_pos + len < seq.size())
        ++len;

    for (size_t i = 0; i < len; ++i) {
        bool is_lower = islower(seq[start_pos + i]);
        ContainCorrected   |= !is_lower;
        ContainUncorrected |=  is_lower;
    }

    if (!ContainUncorrected) // => ContainCorrected == true
        return StatType::on_corrected;
    if (!ContainCorrected)   // => ContainUncorrected == true
        return StatType::on_uncorrected;
    return StatType::on_bound;
}

struct ErrorStatistics {
    Stat mismatch;
    Stat insertion; // to contig
    Stat deletion;  // from contig
};

#define WISHED_COLUMNS Columns::contig_name,\
    Columns::ref_nucls, Columns::contig_nucls, Columns::contig_start_pos

struct FullErrorStatistics {
    ErrorStatistics events;
    ErrorStatistics total_len;

    void AddStatistics(Record<Columns, WISHED_COLUMNS> const & record, std::string const & seq);
};

std::ostream & operator << (std::ostream & out, Stat const & stat) {
    auto total = std::max(std::accumulate(stat.begin(), stat.end(), 0.0), 1.0);
    auto Print = [total, &out] (char const * type, size_t value) {
        out << "       " << type <<": " << value << " (" << (value * 100.0) / total << "%)" << '\n';
    };

    Print("on uncorrected", stat[StatType::on_uncorrected]);
    Print("on corrected  ", stat[StatType::on_corrected]);
    Print("on bound      ", stat[StatType::on_bound]);
    return out;
}

std::ostream & operator << (std::ostream & out, ErrorStatistics const & stat) {
    out << "     mismatches:\n";
    out << stat.mismatch;
    out << "     insertions:\n";
    out << stat.insertion;
    out << "     deletions:\n";
    out << stat.deletion;
    return out;
}

std::ostream & operator << (std::ostream & out, FullErrorStatistics const & stat) {
    out << "   events:\n";
    out << stat.events;
    out << "   total len:\n";
    out << stat.total_len;
    return out;
}

std::string ToUpperStr(std::string && s) {
    for (auto & c : s)
        c = toupper(c);
    return std::move(s);
}

void FullErrorStatistics::AddStatistics(Record<Columns, WISHED_COLUMNS> const & record, std::string const & seq) {
    auto start_pos = record.Get<Columns::contig_start_pos>();
    auto ref_nucls = record.Get<Columns::ref_nucls>();
    auto contig_nucls = record.Get<Columns::contig_nucls>();

    VERIFY(start_pos > 0); // numbering from 1
    --start_pos;

    if (contig_nucls != ".") {
        auto seq_nucls = ToUpperStr(seq.substr(start_pos, contig_nucls.size()));
        VERIFY(seq_nucls == contig_nucls);
    }

    VERIFY(ref_nucls != "." || contig_nucls != ".");
    VERIFY(ref_nucls == "." || contig_nucls == "." || ref_nucls.size() == contig_nucls.size());

    auto len = (contig_nucls != "." ? contig_nucls.size() : 0);
    auto stat_type = GetStatType(seq, start_pos, len);

    if (ref_nucls != "." && contig_nucls != ".") {
        ++events.mismatch[stat_type];
        total_len.mismatch[stat_type] += contig_nucls.size();
    } else if (ref_nucls != ".") {
        ++events.deletion[stat_type];
        total_len.deletion[stat_type] += ref_nucls.size();
    } else {
        ++events.insertion[stat_type];
        total_len.insertion[stat_type] += contig_nucls.size();
    }
}

} // namespace

int main(int argc, char * argv[]) {
    auto contigs = ReadContigs(cfg.contigs);
    auto snps = ReadUsedSnps<WISHED_COLUMNS>(cfg.snps);

    auto index = IndexByName(contigs);

    std::unordered_set<std::string> ignored_contigs;
    std::unordered_map<std::string, FullErrorStatistics> statistics;

    for (auto const & snp : snps) {
        auto const & contig_name = snp.Get<Columns::contig_name>();
        auto id = index.find(contig_name);
        if (id == index.end())
            ignored_contigs.insert(contig_name);

        auto const & seq = contigs[id->second].seq;
        statistics[contig_name].AddStatistics(snp, seq);
    }

    for (auto const & contig : ignored_contigs) 
        std::cout << "The sequence '" + contig +  "' is in file with used snps, but not in file with contigs; ignoring\n";

    for (auto const & stat : statistics) {
        std::cout << '>' << stat.first << '\n';
        std::cout << stat.second;
    }
    return 0;
}

} // namespace error_analyzer
