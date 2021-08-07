#pragma once

#include "common_enums.hpp"
#include "helpers/template_utils.hpp"

#include <array>
#include <cstddef>
#include <numeric>
#include <iosfwd>

namespace error_analyzer {

enum class StatType {
    on_uncorrected = 0,
    on_corrected   = 1,
    on_bound       = 2
};

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

    size_t Sum() const noexcept {
        return std::accumulate(this->begin(), this->end(), size_t(0));
    }
};

using LocalErrorStatType = Stat<3, StatType>;

using CoverageStatistics = Stat<5, RangeType, RangeEndsType>;

struct ErrorStatistics {
    LocalErrorStatType mismatch;
    LocalErrorStatType insertion; // to contig
    LocalErrorStatType deletion;  // from contig
};

struct FullErrorStatistics {
    ErrorStatistics events;
    ErrorStatistics total_len;
    CoverageStatistics cov_stats;
};

std::ostream & operator << (std::ostream & out, CoverageStatistics const & stat);

std::ostream & operator << (std::ostream & out, LocalErrorStatType const & stat);

std::ostream & operator << (std::ostream & out, ErrorStatistics const & stat);

std::ostream & operator << (std::ostream & out, FullErrorStatistics const & stat);

} // namespace error_analyzer
