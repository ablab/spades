#pragma once

#include "common.hpp"
#include <vector>

namespace emphf {

    class bitpair_vector {
    public:
        bitpair_vector(): m_size(0) {}
        bitpair_vector(uint64_t n): m_size(0){resize(n);}
        void resize(uint64_t n);
        size_t size() const;
        size_t mem_size() const;
        uint64_t operator[](uint64_t pos) const;
        void set(uint64_t pos, uint64_t val);
        uint64_t range_nonzeros(uint64_t begin, uint64_t end) const;
        void swap(bitpair_vector& other);
        void save(std::ostream& os) const;
        void load(std::istream& is);
        std::vector<uint64_t> const & data() const;
    protected:
        std::vector<uint64_t> m_bits;
        uint64_t m_size;
    };

}
