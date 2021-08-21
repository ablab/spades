#include "error_analyzer.hpp"
#include "statistics.hpp"
#include "sequence_fragments.hpp"
#include "helpers/common.hpp"
#include "helpers/aligner_output_reader.hpp"
#include "helpers/replace_info_reader.hpp"

#include "utils/filesystem/file_opener.hpp"

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
    auto inp = fs::open_file(file, std::ios_base::in, std::ios_base::badbit);
    Records<SNPSColumns, columns ...> records;
    RecordPusher<SNPSColumns, columns ...> pusher(records, [](...){ return true; });
    std::string line;
    while (GetNextNonemptyLine(inp, line))
        pusher.Push(line);
    return records;
}

#define SNPS_WISHED_COLUMNS SNPSColumns::contig_name, SNPSColumns::ref_name,\
    SNPSColumns::ref_nucls, SNPSColumns::contig_nucls, SNPSColumns::contig_start_pos

#define REPLACE_INFO_WISHED_COLUMNS ReplaceInfoColumns::contig_name, ReplaceInfoColumns::seq_type, \
    ReplaceInfoColumns::from, ReplaceInfoColumns::len

void addCoverageStats(SeqFragments const &fragments, FullErrorStatistics & stats) {
    for (auto const & fragment : fragments)
        stats.cov_stats[fragment.type] += fragment.len;

    if (fragments.empty())
        return;

    stats.cov_stats[RangeEndsType::origin_head] += fragments[0].len;
    if (fragments.size() > 1)
        stats.cov_stats[RangeEndsType::origin_tail] += fragments.back().len;
}

bool IsOnBound(SeqFragments const &fragments, size_t start_pos, size_t end_pos) {
    if (start_pos == end_pos)
        return false;

    for (; start_pos <= end_pos; ++start_pos) {
        if (fragments[start_pos].type == RangeType::origin)
            return true; // because origin fragment always around edges
    }

    return false;
}

void AddBoundStatistics(Record<SNPSColumns, SNPS_WISHED_COLUMNS> const & record, SeqFragments const & seqFragments, FullErrorStatistics & stats) {
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

    if (!IsOnBound(seqFragments, first_and_last_fragments.first, first_and_last_fragments.second))
        return;

    size_t real_len = std::max(contig_nucls.size(), ref_nucls.size());
    ErrorType type;

    if (ref_nucls != "." && contig_nucls != ".")
        type = ErrorType::mismatch;
    else if (ref_nucls != ".")
        type = ErrorType::deletion;
    else
        type = ErrorType::insertion;

    ++stats.events[type][BoundStatType::on_bound];
    stats.total_len[type][BoundStatType::on_bound] += real_len;

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

    auto first_and_last_fragments = seqFragments.GetAllIntersectedFragments(start_pos, len);

    size_t real_len = std::max(contig_nucls.size(), ref_nucls.size());
    ErrorType type;

    if (ref_nucls != "." && contig_nucls != ".")
        type = ErrorType::mismatch;
    else if (ref_nucls != ".")
        type = ErrorType::deletion;
    else
        type = ErrorType::insertion;

    if (first_and_last_fragments.first == first_and_last_fragments.second) {
        auto stat_type = seqFragments[first_and_last_fragments.first].type;
        ++stats.events[type][stat_type];
        stats.total_len[type][stat_type] += real_len;
    } else {
        VERIFY(type != ErrorType::deletion);
        auto end_pos = len + start_pos;
        for (auto i = first_and_last_fragments.first; i <= first_and_last_fragments.second; ++i) {
            auto const & fragment = seqFragments[i];
            auto intersection_len = std::min(fragment.start_pos + fragment.len, end_pos) - start_pos;
            start_pos += intersection_len;
            VERIFY(start_pos <= end_pos);
            ++stats.events[type][fragment.type];
            stats.total_len[type][fragment.type] += intersection_len;
        }
    }

    AddBoundStatistics(record, seqFragments, stats);
}

template <ReplaceInfoColumns ... columns>
std::unordered_map<std::string, SeqFragments> MakeSeqFragments(Records<ReplaceInfoColumns, columns ...> const & records) {
    std::unordered_map<std::string, SeqFragments> res;
    for (auto const & record : records) {
        auto & seqFragments = res[record.template Get<ReplaceInfoColumns::contig_name>()];
        auto len = record.template Get<ReplaceInfoColumns::len>();
        seqFragments.total_len += len;
        auto from = record.template Get<ReplaceInfoColumns::from>();
        auto type = record.template Get<ReplaceInfoColumns::seq_type>();
        seqFragments.push_back({from, len, type});
    }

    for (auto & fragments : res) {
        std::sort(fragments.second.begin(), fragments.second.end());
        if (!fragments.second.isConsequent())
            throw "Replace info for '" + fragments.first + "' is corrupted!";
    }

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

    StatisticsByReference statistics;

    for (auto const & snp : snps) {
        auto const & contig_name = snp.Get<SNPSColumns::contig_name>();
        auto seqFragments = seqFragmentsByName.find(contig_name);
        if (seqFragments == seqFragmentsByName.end())
            continue;

        auto const & ref_name = snp.Get<SNPSColumns::ref_name>();
        AddStatistics(snp, seqFragments->second, statistics[ref_name][contig_name]);
    }

    for (auto & ref_stat : statistics) {
        for (auto & stat : ref_stat.second) {
            auto seqFragments = seqFragmentsByName.find(stat.first);
            addCoverageStats(seqFragments->second, stat.second);
            ref_stat.second.summary_statistics += stat.second;
        }
    }

    for (auto & ref_stat : statistics) {
        std::cout << "----[ " << ref_stat.first << " ]----" << '\n';
        std::cout << "----< summary >----\n";
        std::cout << ref_stat.second.summary_statistics;
        std::cout << "----< end of summary >----\n";
        std::cout << "----{ details }----\n";
        for (auto const & stat : ref_stat.second) {
            std::cout << '>' << stat.first << '\n';
            std::cout << stat.second;
        }
        std::cout << "----{ end of details }----\n";
        std::cout << "----[ end of " << ref_stat.first << " ]----" << '\n';
    }
    return 0;
}

} // namespace error_analyzer
