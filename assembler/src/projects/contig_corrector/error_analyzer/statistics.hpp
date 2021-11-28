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
template<size_t amount, class Element, class ... ValidAccessors>
struct Stat : std::array<Element, amount> {
    using Base = std::array<Element, amount>;

    Stat() : Base({}) {};

    template<class EnumClassType>
    Element & operator[](EnumClassType type) {
        static_assert(traits::Contains<EnumClassType, ValidAccessors ...>());
        return Base::operator[](static_cast<size_t>(type));
    }
    
    template<class EnumClassType>
    Element const & operator[](EnumClassType type) const {
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

using RangeNumberStatType = Stat<3, std::size_t, RangeType>;

using LocalErrorStatType = Stat<3+1, std::size_t, RangeType, BoundStatType>;

using CoverageStatistics = Stat<3+2, std::size_t, RangeType, RangeEndsType>;

using ErrorStatistics = Stat<3, LocalErrorStatType, ErrorType>;

struct FullErrorStatistics {
    ErrorStatistics events;
    ErrorStatistics total_len;
    CoverageStatistics cov_stats;
    RangeNumberStatType range_num_stat;

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
