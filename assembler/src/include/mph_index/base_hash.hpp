#pragma once

#include <cstdint>
#include <tuple>
#include <algorithm>
#include <cstring>
#include "common.hpp"


namespace emphf {

    inline uint64_t unaligned_load64(uint8_t const* from)
    {
        uint64_t tmp;
        memcpy(reinterpret_cast<char*>(&tmp), from, 8);
        // XXX(ot): reverse bytes in big-endian architectures
        return tmp;
    }


    struct jenkins64_hasher {

        typedef uint64_t seed_t;
        typedef uint64_t hash_t;
        typedef std::tuple<hash_t, hash_t, hash_t> hash_triple_t;

        jenkins64_hasher()
        {}

        jenkins64_hasher(uint64_t seed)
            : m_seed(seed)
        {}

        template <typename Rng>
        static jenkins64_hasher generate(Rng& rng)
        {
            return jenkins64_hasher(rng());
        }

        // Adapted from http://www.burtleburtle.net/bob/c/lookup8.c
        hash_triple_t operator()(byte_range_t s) const
        {
            using std::get;
            hash_triple_t h(m_seed, m_seed, 0x9e3779b97f4a7c13ULL);

            size_t len = (size_t)(s.second - s.first);
            uint8_t const* cur = s.first;
            uint8_t const* end = s.second;

            while (end - cur >= 24) {
                get<0>(h) += unaligned_load64(cur);
                cur += 8;
                get<1>(h) += unaligned_load64(cur);
                cur += 8;
                get<2>(h) += unaligned_load64(cur);
                cur += 8;

                mix(h);
            }

            get<2>(h) += len;

            switch (end - cur) {
            case 23: get<2>(h) += (uint64_t(cur[22]) << 56);
            case 22: get<2>(h) += (uint64_t(cur[21]) << 48);
            case 21: get<2>(h) += (uint64_t(cur[20]) << 40);
            case 20: get<2>(h) += (uint64_t(cur[19]) << 32);
            case 19: get<2>(h) += (uint64_t(cur[18]) << 24);
            case 18: get<2>(h) += (uint64_t(cur[17]) << 16);
            case 17: get<2>(h) += (uint64_t(cur[16]) << 8);
                // the first byte of c is reserved for the length
            case 16: get<1>(h) += (uint64_t(cur[15]) << 56);
            case 15: get<1>(h) += (uint64_t(cur[14]) << 48);
            case 14: get<1>(h) += (uint64_t(cur[13]) << 40);
            case 13: get<1>(h) += (uint64_t(cur[12]) << 32);
            case 12: get<1>(h) += (uint64_t(cur[11]) << 24);
            case 11: get<1>(h) += (uint64_t(cur[10]) << 16);
            case 10: get<1>(h) += (uint64_t(cur[ 9]) << 8);
            case  9: get<1>(h) += (uint64_t(cur[ 8]));
            case  8: get<0>(h) += (uint64_t(cur[ 7]) << 56);
            case  7: get<0>(h) += (uint64_t(cur[ 6]) << 48);
            case  6: get<0>(h) += (uint64_t(cur[ 5]) << 40);
            case  5: get<0>(h) += (uint64_t(cur[ 4]) << 32);
            case  4: get<0>(h) += (uint64_t(cur[ 3]) << 24);
            case  3: get<0>(h) += (uint64_t(cur[ 2]) << 16);
            case  2: get<0>(h) += (uint64_t(cur[ 1]) << 8);
            case  1: get<0>(h) += (uint64_t(cur[ 0]));
            case 0: break; // nothing to add
            default: assert(false);
            }

            mix(h);

            return h;
        }

        // rehash a hash triple
        hash_triple_t operator()(hash_triple_t h) const
        {
            std::get<0>(h) += m_seed;
            std::get<1>(h) += m_seed;
            std::get<2>(h) += 0x9e3779b97f4a7c13ULL;

            mix(h);

            return h;
        }

        void swap(jenkins64_hasher& other)
        {
            std::swap(m_seed, other.m_seed);
        }

        void save(std::ostream& os) const
        {
            os.write(reinterpret_cast<char const*>(&m_seed), sizeof(m_seed));
        }

        void load(std::istream& is)
        {
            is.read(reinterpret_cast<char*>(&m_seed), sizeof(m_seed));
        }

        seed_t seed() const
        {
            return m_seed;
        }

    protected:

        static void mix(hash_triple_t& h)
        {
            uint64_t& a = std::get<0>(h);
            uint64_t& b = std::get<1>(h);
            uint64_t& c = std::get<2>(h);

            a -= b; a -= c; a ^= (c >> 43);
            b -= c; b -= a; b ^= (a << 9);
            c -= a; c -= b; c ^= (b >> 8);
            a -= b; a -= c; a ^= (c >> 38);
            b -= c; b -= a; b ^= (a << 23);
            c -= a; c -= b; c ^= (b >> 5);
            a -= b; a -= c; a ^= (c >> 35);
            b -= c; b -= a; b ^= (a << 49);
            c -= a; c -= b; c ^= (b >> 11);
            a -= b; a -= c; a ^= (c >> 12);
            b -= c; b -= a; b ^= (a << 18);
            c -= a; c -= b; c ^= (b >> 22);
        }

        seed_t m_seed;
    };


    // This is basically a wrapper to jenkins64_hasher that uses a
    // 32-bit seed and returns 32-bit hashes by truncation
    struct jenkins32_hasher {

        typedef uint32_t seed_t;
        typedef uint32_t hash_t;
        typedef std::tuple<hash_t, hash_t, hash_t> hash_triple_t;

        jenkins32_hasher()
        {}

        jenkins32_hasher(uint32_t seed)
            : m_seed(seed)
        {}

        template <typename Rng>
        static jenkins32_hasher generate(Rng& rng)
        {
            return jenkins32_hasher((uint32_t)rng());
        }

        hash_triple_t operator()(byte_range_t s) const
        {
            using std::get;
            auto h64 = jenkins64_hasher(seed64())(s);
            return hash_triple_t((uint32_t)get<0>(h64),
                                 (uint32_t)get<1>(h64),
                                 (uint32_t)get<2>(h64));
        }

        hash_triple_t operator()(hash_triple_t h) const
        {
            using std::get;
            auto h64 = jenkins64_hasher::hash_triple_t(get<0>(h),
                                                       get<1>(h),
                                                       get<2>(h));
            h64 = jenkins64_hasher(seed64())(h64);
            return hash_triple_t((uint32_t)get<0>(h64),
                                 (uint32_t)get<1>(h64),
                                 (uint32_t)get<2>(h64));
        }

        void swap(jenkins32_hasher& other)
        {
            std::swap(m_seed, other.m_seed);
        }

        void save(std::ostream& os) const
        {
            os.write(reinterpret_cast<char const*>(&m_seed), sizeof(m_seed));
        }

        void load(std::istream& is)
        {
            is.read(reinterpret_cast<char*>(&m_seed), sizeof(m_seed));
        }

        seed_t seed() const
        {
            return m_seed;
        }

    protected:

        uint64_t seed64() const
        {
            return (uint64_t(m_seed) << 32) | m_seed;
        }

        seed_t m_seed;

    };

}
