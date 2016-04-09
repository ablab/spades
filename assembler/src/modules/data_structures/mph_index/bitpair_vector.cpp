//
// Created by anton on 3/22/16.
//

#include "bitpair_vector.hpp"

#include <iostream>

void emphf::bitpair_vector::resize(uint64_t n) {
    // can only grow, for now
    assert(n >= size());
    m_size = n;
    m_bits.resize((m_size + 31) / 32);
}

size_t emphf::bitpair_vector::size() const {
    return m_size;
}

size_t emphf::bitpair_vector::mem_size() const {
    return m_bits.size() * sizeof(m_bits[0]);
}

uint64_t emphf::bitpair_vector::operator[](uint64_t pos) const {
    return (m_bits[pos / 32] >> ((pos % 32) * 2)) % 4;
}

void emphf::bitpair_vector::set(uint64_t pos, uint64_t val) {
    assert(val < 4);
    uint64_t word_pos = pos / 32;
    uint64_t word_offset = (pos % 32) * 2;
    m_bits[word_pos] &= ~(3ULL << word_offset);
    m_bits[word_pos] |= val << word_offset;
}

uint64_t emphf::bitpair_vector::range_nonzeros(uint64_t begin, uint64_t end) const {
    assert(begin <= end);
    assert(end <= size());

    uint64_t word_begin = begin / 32;
    uint64_t offset_begin = (begin % 32) * 2;
    uint64_t word_end = end / 32;
    uint64_t offset_end = (end % 32) * 2;
    uint64_t r = 0;

    uint64_t word = (m_bits[word_begin] >> offset_begin) << offset_begin;
    for (uint64_t w = word_begin; w < word_end; ++w) {
        r += nonzero_pairs(word);
        word = m_bits[w + 1];
    }

    uint64_t mask = (uint64_t(1) << offset_end) - 1;
    r += nonzero_pairs(word & mask);

    return r;
}

void emphf::bitpair_vector::swap(bitpair_vector& other) {
            std::swap(m_size, other.m_size);
            m_bits.swap(other.m_bits);
        }


void emphf::bitpair_vector::save(std::ostream& os) const {
            os.write(reinterpret_cast<char const*>(&m_size), sizeof(m_size));
            os.write(reinterpret_cast<char const*>(m_bits.data()), (std::streamsize)(sizeof(m_bits[0]) * m_bits.size()));
        }

void emphf::bitpair_vector::load(std::istream& is) {
            is.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
            m_bits.resize((m_size + 31) / 32);
            is.read(reinterpret_cast<char*>(m_bits.data()), (std::streamsize)(sizeof(m_bits[0]) * m_bits.size()));
        }

std::vector<uint64_t> const &emphf::bitpair_vector::data() const {
    return m_bits;
}
