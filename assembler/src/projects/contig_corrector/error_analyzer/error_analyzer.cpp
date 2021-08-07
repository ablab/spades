#include "error_analyzer.hpp"
#include "statistics.hpp"
#include "sequence_fragments.hpp"
#include "helpers/common.hpp"
#include "helpers/aligner_output_reader.hpp"

#include <iostream>
#include <unordered_map>
#include <unordered_set>

namespace helpers {

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

} // namespace helpers

namespace error_analyzer {

using namespace helpers;

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

#define SNPS_WISHED_COLUMNS SNPSColumns::contig_name,\
    SNPSColumns::ref_nucls, SNPSColumns::contig_nucls, SNPSColumns::contig_start_pos

#define REPLACE_INFO_WISHED_COLUMNS ReplaceInfoColumns::contig_name, ReplaceInfoColumns::seq_type, \
    ReplaceInfoColumns::from, ReplaceInfoColumns::len

void addCoverageStats(SeqFragments const &fragments, FullErrorStatistics & stats) {
    for (size_t i = 0; i + 1 < fragments.size(); ++i)
        stats.cov_stats[fragments[i].type] += fragments[i+1].start_pos - fragments[i].start_pos;

    if (fragments.empty())
        return;

    stats.cov_stats[fragments.back().type] += fragments.total_len - fragments.back().start_pos;

    VERIFY_MSG(stats.cov_stats.Sum() == fragments.total_len, "Replace info was corrupted");

    stats.cov_stats[RangeEndsType::origin_tail] += fragments.total_len - fragments.back().start_pos;

    if (fragments.size() > 1)
        stats.cov_stats[RangeEndsType::origin_head] += fragments[1].start_pos - fragments[0].start_pos;
}

StatType GetStatType(SeqFragments const &fragments, size_t start_pos, size_t end_pos) {
    if (start_pos == end_pos) {
        if (fragments[start_pos].type == RangeType::origin)
            return StatType::on_uncorrected;
        return StatType::on_corrected;
    }

    for (; start_pos <= end_pos; ++start_pos) {
        if (fragments[start_pos].type == RangeType::origin)
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
        seqFragments.emplace_back(from, type);
    }

    for (auto & fragments : res)
        std::sort(fragments.second.begin(), fragments.second.end());

    return res;
}

} // namespace

int main() {
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

    for (auto & stat : statistics) {
        auto seqFragments = seqFragmentsByName.find(stat.first);
        addCoverageStats(seqFragments->second, stat.second);
    }

    for (auto const & stat : statistics) {
        std::cout << '>' << stat.first << '\n';
        std::cout << stat.second;
    }
    return 0;
}

} // namespace error_analyzer
