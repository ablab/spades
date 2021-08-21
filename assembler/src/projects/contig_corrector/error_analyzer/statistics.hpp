#pragma once

#include "common_enums.hpp"
#include "helpers/template_utils.hpp"

#include <array>
#include <unordered_map>
#include <string>
#include <cstddef>
#include <iosfwd>
#include <cassert>

namespace error_analyzer {
template<size_t amount, class ... ValidAccessors>
struct Stat : std::array<size_t, amount> {
    using Base = std::array<size_t, amount>;

    Stat() : Base({0}) {};

    template<class EnumClassType>
    size_t & operator[](EnumClassType type) {
        static_assert(traits::Contains<EnumClassType, ValidAccessors ...>());
        return Base::operator[](static_cast<size_t>(type));
    }
    
    template<class EnumClassType>
    size_t const & operator[](EnumClassType type) const {
        static_assert(traits::Contains<EnumClassType, ValidAccessors ...>());
        return Base::operator[](static_cast<size_t>(type));
    }

    void operator += (Stat const & other) {
        auto& this_base = *static_cast<Base*>(this);
        auto& other_base = *static_cast<Base const *>(&other);
        for (size_t i = 0; i < amount; ++i)
            this_base[i] += other_base[i];
    }
};

using LocalErrorStatType = Stat<3+1, RangeType, BoundStatType>;

using CoverageStatistics = Stat<3+2, RangeType, RangeEndsType>;


struct ErrorStatistics {
    LocalErrorStatType mismatch;
    LocalErrorStatType insertion; // to contig
    LocalErrorStatType deletion;  // from contig

    LocalErrorStatType& operator[] (ErrorType type) {
        switch (type) {
            case ErrorType::mismatch : return mismatch; break;
            case ErrorType::insertion: return insertion; break;
            case ErrorType::deletion : return deletion; break;
            default: assert(false && "unreachable");
        }
    }

    void operator +=(ErrorStatistics const & other) {
        mismatch += other.mismatch;
        insertion += other.insertion;
        deletion += other.deletion;
    }
};

struct FullErrorStatistics {
    ErrorStatistics events;
    ErrorStatistics total_len;
    CoverageStatistics cov_stats;

    void operator +=(FullErrorStatistics const & other);
};

std::ostream & operator << (std::ostream & out, CoverageStatistics const & stat);

std::ostream & operator << (std::ostream & out, LocalErrorStatType const & stat);

std::ostream & operator << (std::ostream & out, ErrorStatistics const & stat);

std::ostream & operator << (std::ostream & out, FullErrorStatistics const & stat);

struct StatisticsByContig : std::unordered_map<std::string, FullErrorStatistics> {
    FullErrorStatistics summary_statistics;
};

using StatisticsByReference = std::unordered_map<std::string, StatisticsByContig>;

} // namespace error_analyzer
