#pragma once

#include <cstdint>

#include "common.hpp"
#include "bitpair_vector.hpp"

namespace emphf {

    class ranked_bitpair_vector {
    public:

        ranked_bitpair_vector()
        {}

        void build(bitpair_vector&& bv)
        {
            m_bv.swap(bv);

            uint64_t cur_rank = 0;
            auto const& words = m_bv.data();
            for (size_t i = 0; i < words.size(); ++i) {
                if (((i * 32) % pairs_per_block) == 0) {
                    m_block_ranks.push_back(cur_rank);
                }
                cur_rank += nonzero_pairs(words[i]);
            }
        }

        size_t size() const
        {
            return m_bv.size();
        }

        size_t mem_size() const {
            return m_bv.mem_size() + m_block_ranks.size() * sizeof(m_block_ranks[0]);
        }

        uint64_t operator[](uint64_t pos) const
        {
            return m_bv[pos];
        }

        uint64_t rank(uint64_t pos) const
        {
            uint64_t word_idx = pos / 32;
            uint64_t word_offset = pos % 32;
            uint64_t block = pos / pairs_per_block;
            uint64_t r = m_block_ranks[block];

            for (uint64_t w = block * pairs_per_block / 32; w < word_idx; ++w) {
                r += nonzero_pairs(m_bv.data()[w]);
            }

            uint64_t mask = (uint64_t(1) << (word_offset * 2)) - 1;
            r += nonzero_pairs(m_bv.data()[word_idx] & mask);

            return r;
        }

        void swap(ranked_bitpair_vector& other)
        {
            m_bv.swap(other.m_bv);
            m_block_ranks.swap(other.m_block_ranks);
        }

        void save(std::ostream& os) const
        {
            m_bv.save(os);
            assert(m_block_ranks.size() ==
                   (m_bv.size() + pairs_per_block - 1) / pairs_per_block);
            os.write(reinterpret_cast<char const*>(m_block_ranks.data()),
                     (std::streamsize)(sizeof(m_block_ranks[0]) * m_block_ranks.size()));
        }

        void load(std::istream& is)
        {
            m_bv.load(is);
            m_block_ranks.resize((m_bv.size() + pairs_per_block - 1) / pairs_per_block);
            is.read(reinterpret_cast<char*>(m_block_ranks.data()),
                    (std::streamsize)(sizeof(m_block_ranks[0]) * m_block_ranks.size()));
        }

    protected:

        static const uint64_t pairs_per_block = 512;
        bitpair_vector m_bv;
        std::vector<uint64_t> m_block_ranks;
    };

}
