#include "statistics.hpp"

#include <fstream>

namespace error_analyzer {

std::ostream & operator << (std::ostream & out, CoverageStatistics const & stat) {
    auto total = std::max<size_t>(stat[RangeType::origin] + stat[RangeType::edge] + stat[RangeType::path], 1);
    auto Print = [total, &out] (char const * type, size_t value) {
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
    auto total = std::max((double)stat.Sum(), 1.0);
    auto Print = [total, &out] (char const * type, size_t value) {
        out << "       " << type <<": " << value << " (" << (static_cast<double>(value) * 100.0) / total << "%)" << '\n';
    };

    Print("on uncorrected", stat[StatType::on_uncorrected]);
    Print("on corrected  ", stat[StatType::on_corrected]);
    Print("on bound      ", stat[StatType::on_bound]);
    return out;
}

std::ostream & operator << (std::ostream & out, ErrorStatistics const & stat) {
    out << "     mismatches: " << stat.mismatch.Sum() << " \n";
    out << stat.mismatch;
    out << "     insertions: " << stat.insertion.Sum() << " \n";
    out << stat.insertion;
    out << "     deletions: " << stat.deletion.Sum() << " \n";
    out << stat.deletion;
    return out;
}

std::ostream & operator << (std::ostream & out, FullErrorStatistics const & stat) {
    out << "   coverage fragment length:\n";
    out << stat.cov_stats;
    out << "   events:\n";
    out << stat.events;
    out << "   total len:\n";
    out << stat.total_len;
    return out;
}


} // namespace error_analyzer
