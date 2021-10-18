#include "statistics.hpp"

#include <fstream>

namespace error_analyzer {

template<size_t s, class ... Args>
size_t RangeTypeSum(Stat<s, Args...> const & stat) {
    return stat[RangeType::origin] + stat[RangeType::edge] + stat[RangeType::path];
}

void FullErrorStatistics::operator +=(FullErrorStatistics const & other) {
    events += other.events;
    total_len += other.total_len;
    cov_stats += other.cov_stats;
    range_num_stat += other.range_num_stat;
}

std::ostream & operator << (std::ostream & out, RangeNumberStatType const & stat) {
    size_t total = RangeTypeSum(stat);
    if (!total)
        return out;
    auto Print = [total, &out] (char const * type, size_t value) {
        if (value)
            out << "     " << type <<": " << value << " (" << (static_cast<double>(value) * 100.0) / static_cast<double>(total) << "%)" << '\n';
    };

    Print("origin", stat[RangeType::origin]);
    Print("edges ", stat[RangeType::edge]);
    Print("paths ", stat[RangeType::path]);
    return out;
}

std::ostream & operator << (std::ostream & out, CoverageStatistics const & stat) {
    size_t total = RangeTypeSum(stat);
    if (!total)
        return out;
    auto Print = [total, &out] (char const * type, size_t value) {
        if (value)
            out << "     " << type <<": " << value << " (" << (static_cast<double>(value) * 100.0) / static_cast<double>(total) << "%)" << '\n';
    };

    Print("uncorrected", stat[RangeType::origin]);
    Print("  - head ", stat[RangeEndsType::origin_head]);
    Print("  - tail ", stat[RangeEndsType::origin_tail]);
    Print("corrected  ", stat[RangeType::edge] + stat[RangeType::path]);
    Print("  - edges", stat[RangeType::edge]);
    Print("  - paths", stat[RangeType::path]);
    return out;
}

std::ostream & operator << (std::ostream & out, LocalErrorStatType const & stat) {
    size_t total = RangeTypeSum(stat);
    if (!total)
        return out;
    auto Print = [total, &out] (char const * type, size_t value) {
        if (value)
            out << "       " << type <<": " << value << " (" << (static_cast<double>(value) * 100.0) / static_cast<double>(total) << "%)" << '\n';
    };

    Print("on uncorrected", stat[RangeType::origin]);
    Print("on corrected  ", stat[RangeType::edge] + stat[RangeType::path]);
    Print("  - edges ", stat[RangeType::edge]);
    Print("  - paths ", stat[RangeType::path]);
    Print("on bounds     ", stat[BoundStatType::on_bound]);
    return out;
}

std::ostream & operator << (std::ostream & out, ErrorStatistics const & stat) {
    out << "     mismatches: " << RangeTypeSum(stat.mismatch) << " \n";
    out << stat.mismatch;
    out << "     insertions: " << RangeTypeSum(stat.insertion) << " \n";
    out << stat.insertion;
    out << "     deletions: " << RangeTypeSum(stat.deletion) << " \n";
    out << stat.deletion;
    return out;
}

std::ostream & operator << (std::ostream & out, FullErrorStatistics const & stat) {
    out << "   fragment number:\n";
    out << stat.range_num_stat;
    out << "   coverage fragment length:\n";
    out << stat.cov_stats;
    out << "   events:\n";
    out << stat.events;
    out << "   total len:\n";
    out << stat.total_len;
    return out;
}


} // namespace error_analyzer
