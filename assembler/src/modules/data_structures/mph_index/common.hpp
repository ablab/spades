#pragma once

#include <cstdint>
#include <iterator>
#include <memory>
#include <cassert>

#include "emphf_config.hpp"

namespace emphf {

    template <typename Iterator>
    struct iter_range
    {
        iter_range(Iterator b, Iterator e)
            : m_begin(b)
            , m_end(e)
        {}

        Iterator begin() const
        { return m_begin; }

        Iterator end() const
        { return m_end; }

        Iterator m_begin, m_end;
    };

    typedef std::pair<uint8_t const*, uint8_t const*> byte_range_t;

    struct identity_adaptor
    {
        byte_range_t operator()(byte_range_t s) const
        {
            return s;
        }
    };

    template <typename Iterator>
    iter_range<Iterator> range(Iterator begin, Iterator end)
    {
        return iter_range<Iterator>(begin, end);
    }

    inline uint64_t nonzero_pairs(uint64_t x)
    {
        static const uint64_t ones_step_4  = 0x1111111111111111ULL;
        x = (x | (x >> 1)) & (0x5 * ones_step_4);

#if EMPHF_USE_POPCOUNT
        return (uint64_t)__builtin_popcountll(x);
#else
        static const uint64_t ones_step_8  = 0x0101010101010101ULL;
        x = (x & 3 * ones_step_4) + ((x >> 2) & 3 * ones_step_4);
        x = (x + (x >> 4)) & 0x0f * ones_step_8;
        return (x * ones_step_8) >> 56;
#endif
    }

    inline uint64_t msb(uint64_t x)
    {
        assert(x);
        return 63 - __builtin_clzll(x);
    }

}
