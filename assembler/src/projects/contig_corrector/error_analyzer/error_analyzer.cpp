#include "error_analyzer.hpp"
#include "helpers/common.hpp"
#include "helpers/aligner_output_reader.hpp"

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <numeric>

enum class SNPSColumns : size_t {
    ref_name           = 0,
    contig_name        = 1,
    ref_start_pos      = 2,
    ref_nucls          = 3,
    contig_nucls       = 4,
    contig_start_pos   = 5,
    TOTAL_COLUMNS_SIZE = 6
};

template<SNPSColumns el>
struct type_getter<SNPSColumns, el> {
    using type = std::conditional_t<
                        el == SNPSColumns::ref_start_pos ||
                        el == SNPSColumns::contig_start_pos,
                        unsigned long long,
                        std::string>;
};

enum class ReplaceInfoColumns : size_t {
    contig_name        = 0,
    seq_type           = 1,
    from               = 2,
    len                = 3,
    TOTAL_COLUMNS_SIZE = 4
};

template<ReplaceInfoColumns el>
struct type_getter<ReplaceInfoColumns, el> {
    using type = std::conditional_t<
                        el == ReplaceInfoColumns::from ||
                        el == ReplaceInfoColumns::len,
                        unsigned long long,
                        std::string>;
};

namespace error_analyzer {

struct gcfg {
    std::string replace_info;
    std::string snps;
} cfg;


clipp::group GetCLI() {
  using namespace clipp;

  auto cli = (
      cfg.replace_info << value("replace info dump file"),
      cfg.snps << value("used_snps.gz")
  );

  return cli;
}

namespace {

template<SNPSColumns ... columns>
Records<SNPSColumns, columns ...> ReadUsedSnps(std::string const & file) {
    std::ifstream inp(file);
    if (!inp.is_open())
        throw "Cannot open '" + file + "' file";
    Records<SNPSColumns, columns ...> records;
    RecordPusher<SNPSColumns, columns ...> pusher(records, [](...){ return true; });
    std::string line;
    while (GetNextNonemptyLine(inp, line))
        pusher.Push(line);
    return records;
}

template<ReplaceInfoColumns ... columns>
Records<ReplaceInfoColumns, columns ...> ReadReplaceInfoDump(std::string const & file, FilterType<ReplaceInfoColumns, columns ...> filter) {
    std::ifstream inp(file);
    if (!inp.is_open())
        throw "Cannot open '" + file + "' file";
    Records<ReplaceInfoColumns, columns ...> records;
    RecordPusher<ReplaceInfoColumns, columns ...> pusher(records, filter);
    std::string line;
    GetNextNonemptyLine(inp, line); // skipping header
    while (GetNextNonemptyLine(inp, line))
        pusher.Push(line);
    return records;
}

enum class RangeType : char {
    origin = 0,
    edge   = 1,
    path   = 2
};

enum class StatType {
    on_uncorrected = 0,
    on_corrected   = 1,
    on_bound       = 2
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

struct ErrorStatistics {
    Stat mismatch;
    Stat insertion; // to contig
    Stat deletion;  // from contig
};

#define SNPS_WISHED_COLUMNS SNPSColumns::contig_name,\
    SNPSColumns::ref_nucls, SNPSColumns::contig_nucls, SNPSColumns::contig_start_pos

#define REPLACE_INFO_WISHED_COLUMNS ReplaceInfoColumns::contig_name, ReplaceInfoColumns::seq_type, \
    ReplaceInfoColumns::from, ReplaceInfoColumns::len

struct FullErrorStatistics {
    ErrorStatistics events;
    ErrorStatistics total_len;
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

template<class Columns, Columns field_name, Columns ... columns>
std::unordered_set<type_getter_t<Columns, field_name>> CollectUniqueFieldValues(Records<Columns, columns ...> const & snps){
    std::unordered_set<type_getter_t<Columns, field_name>> res;
    for (auto const & snp : snps)
        res.insert(snp.template Get<field_name>());
    return res;
}

struct SeqFragment {
    unsigned long long start_pos : 62;
    RangeType type : 2;

    SeqFragment(unsigned long long from, RangeType type)
        : start_pos(from)
        , type(type)
    {}

    bool operator < (SeqFragment const & other) const noexcept {
        return start_pos < other.start_pos;
    }
};

struct SeqFragments {
    std::vector<SeqFragment> fragments;
    size_t total_len = 0;

    size_t FindFragmentIndex(unsigned long long pos) const noexcept {
        size_t l = 0, r = fragments.size();
        while (l + 1 < r) {
            auto mid = (l + r) / 2;
            if (pos < fragments[mid].start_pos)
                r = mid;
            else
                l = mid;
        }
        return l;
    }

    std::pair<size_t, size_t> GetAllIntersectedFragments(unsigned long long start_pos, size_t len) const noexcept {
        auto first_fragment = FindFragmentIndex(start_pos);
        auto end_pos = start_pos + len;
        VERIFY(end_pos <= total_len);
        auto end_fragment = first_fragment + 1;
        while (end_fragment < fragments.size() && fragments[end_fragment].start_pos < end_pos)
            ++end_fragment;
        return {first_fragment, end_fragment - 1};
    }
};

StatType GetStatType(SeqFragments const &fragments, size_t start_pos, size_t end_pos) {
    if (start_pos == end_pos) {
        if (fragments.fragments[start_pos].type == RangeType::origin)
            return StatType::on_uncorrected;
        return StatType::on_corrected;
    }

    for (; start_pos <= end_pos; ++start_pos) {
        if (fragments.fragments[start_pos].type == RangeType::origin)
            return StatType::on_bound; // because origin fragment always around edges
    }

    return StatType::on_corrected;
}

void AddStatistics(Record<SNPSColumns, SNPS_WISHED_COLUMNS> const & record, SeqFragments const & seqFragments, FullErrorStatistics & stats) {
    auto start_pos = record.Get<SNPSColumns::contig_start_pos>();
    auto ref_nucls = record.Get<SNPSColumns::ref_nucls>();
    auto contig_nucls = record.Get<SNPSColumns::contig_nucls>();

    VERIFY(start_pos > 0); // numbering from 1
    --start_pos;

    VERIFY(ref_nucls != "." || contig_nucls != ".");
    VERIFY(ref_nucls == "." || contig_nucls == "." || ref_nucls.size() == contig_nucls.size());

    auto len = (contig_nucls != "." ? contig_nucls.size() : 0);

    // extend checking range to left by 1 if possible
    if (start_pos > 0) {
        --start_pos;
        ++len;
    }

    // extend checking range to right by 1 if possible
    if (start_pos + len < seqFragments.total_len)
        ++len;

    auto first_and_last_fragments = seqFragments.GetAllIntersectedFragments(start_pos, len);

    auto stat_type = GetStatType(seqFragments, first_and_last_fragments.first, first_and_last_fragments.second);

    if (ref_nucls != "." && contig_nucls != ".") {
        ++stats.events.mismatch[stat_type];
        stats.total_len.mismatch[stat_type] += contig_nucls.size();
    } else if (ref_nucls != ".") {
        ++stats.events.deletion[stat_type];
        stats.total_len.deletion[stat_type] += ref_nucls.size();
    } else {
        ++stats.events.insertion[stat_type];
        stats.total_len.insertion[stat_type] += contig_nucls.size();
    }
}

RangeType RangeTypeFromStr(std::string const & s) {
    if (s == "origin")
        return RangeType::origin;
    if (s == "edge")
        return RangeType::edge;
    if (s == "path")
        return RangeType::path;
    throw std::invalid_argument("Unknown type value '" + s + "' in replace info; suppoted values: 'origin'/'edge'/'path'");
}

template <ReplaceInfoColumns ... columns>
std::unordered_map<std::string, SeqFragments> MakeSeqFragments(Records<ReplaceInfoColumns, columns ...> const & records) {
    std::unordered_map<std::string, SeqFragments> res;
    for (auto const & record : records) {
        auto & seqFragments = res[record.template Get<ReplaceInfoColumns::contig_name>()];
        seqFragments.total_len += record.template Get<ReplaceInfoColumns::len>();
        auto from = record.template Get<ReplaceInfoColumns::from>();
        auto type = RangeTypeFromStr(record.template Get<ReplaceInfoColumns::seq_type>());
        seqFragments.fragments.emplace_back(from, type);
    }

    for (auto & fragments : res)
        std::sort(fragments.second.fragments.begin(), fragments.second.fragments.end());

    return res;
}

} // namespace

int main(int argc, char * argv[]) {
    auto snps = ReadUsedSnps<SNPS_WISHED_COLUMNS>(cfg.snps);
    auto snps_contig_names = CollectUniqueFieldValues<SNPSColumns, SNPSColumns::contig_name>(snps);

    auto DropByMissedContigName = [&snps_contig_names](auto const & record) {
            return snps_contig_names.count(record.template Get<ReplaceInfoColumns::contig_name>()) > 0;
        };
    auto replace_info = ReadReplaceInfoDump<REPLACE_INFO_WISHED_COLUMNS>(cfg.replace_info, {DropByMissedContigName});

    auto seqFragmentsByName = MakeSeqFragments(replace_info);

    std::unordered_set<std::string> ignored_contigs;
    auto does_not_contain = [&seqFragmentsByName](auto const & s) { return !seqFragmentsByName.count(s); };
    std::copy_if(snps_contig_names.begin(), snps_contig_names.end(), std::inserter(ignored_contigs, ignored_contigs.end()), does_not_contain);
    for (auto const & contig : ignored_contigs) 
        std::cout << "The sequence '" + contig +  "' is in file with used snps, but not in the replace info dump; ignoring\n";

    std::unordered_map<std::string, FullErrorStatistics> statistics;

    for (auto const & snp : snps) {
        auto const & contig_name = snp.Get<SNPSColumns::contig_name>();
        auto seqFragments = seqFragmentsByName.find(contig_name);
        if (seqFragments == seqFragmentsByName.end())
            continue;

        AddStatistics(snp, seqFragments->second, statistics[contig_name]);
    }

    for (auto const & stat : statistics) {
        std::cout << '>' << stat.first << '\n';
        std::cout << stat.second;
    }
    return 0;
}

} // namespace error_analyzer
