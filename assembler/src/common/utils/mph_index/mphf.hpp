#pragma once

#include <random>

#include "bitpair_vector.hpp"
#include "ranked_bitpair_vector.hpp"

#include "utils/logger/logger.hpp"

namespace emphf {

    template <typename BaseHasher>
    class mphf {
    public:
        mphf()
        {}

        template <typename HypergraphSorter, typename Range, typename Adaptor>
        mphf(HypergraphSorter& sorter, size_t n,
             Range const& input_range, Adaptor adaptor,
             double gamma = 1.23)
            : m_n(n)
            , m_hash_domain(std::max((size_t(std::ceil(double(m_n) * gamma)) + 2) / 3, size_t(2)))
        {
            typedef typename HypergraphSorter::node_t node_t;
            typedef typename HypergraphSorter::hyperedge hyperedge;
            typedef decltype(*std::begin(input_range)) value_type;

            size_t nodes_domain = m_hash_domain * 3;

            if (nodes_domain >= std::numeric_limits<node_t>::max()) {
                throw std::invalid_argument("Too many nodes for node_t");
            }

            auto edge_gen = [&](value_type s) {
                using std::get;
                auto hashes = m_hasher(adaptor(s));
                return hyperedge((node_t)(get<0>(hashes) % m_hash_domain),
                                 (node_t)(m_hash_domain +
                                          (get<1>(hashes) % m_hash_domain)),
                                 (node_t)(2 * m_hash_domain +
                                          (get<2>(hashes) % m_hash_domain)));
            };

            std::mt19937_64 rng(37); // deterministic seed

            for (size_t trial = 0; ; ++trial) {
                //logger() << "Hypergraph generation: trial " << trial << std::endl;

                m_hasher = BaseHasher::generate(rng);
                if (sorter.try_generate_and_sort(input_range, edge_gen,
                                                 m_n, m_hash_domain)) break;
            }

            auto peeling_order = sorter.get_peeling_order();
            bitpair_vector bv(nodes_domain);

            //logger() << "Assigning values" << std::endl;

            for (auto edge = peeling_order.first;
                 edge != peeling_order.second;
                 ++edge) {

                uint64_t target = orientation(*edge);
                uint64_t assigned = bv[edge->v1] + bv[edge->v2];

                // "assigned values" must be nonzeros to be ranked, so
                // if the result is 0 we assign 3
                bv.set(edge->v0, ((target - assigned + 9) % 3) ?: 3);
            }

            m_bv.build(std::move(bv));
        }

        uint64_t size() const
        {
            return m_n;
        }

        size_t mem_size() const {
            return m_bv.mem_size();
        }

        BaseHasher const& base_hasher() const
        {
            return m_hasher;
        }

        template <typename T, typename Adaptor>
        uint64_t lookup(const T &val, Adaptor adaptor)
        {
            using std::get;
            auto hashes = m_hasher(adaptor(val));
            uint64_t nodes[3] = {get<0>(hashes) % m_hash_domain,
                                 m_hash_domain + (get<1>(hashes) % m_hash_domain),
                                 2 * m_hash_domain + (get<2>(hashes) % m_hash_domain)};

            uint64_t hidx = (m_bv[nodes[0]] + m_bv[nodes[1]] + m_bv[nodes[2]]) % 3;
            return m_bv.rank(nodes[hidx]);
        }

        void swap(mphf& other)
        {
            std::swap(m_n, other.m_n);
            std::swap(m_hash_domain, other.m_hash_domain);
            m_hasher.swap(other.m_hasher);
            m_bv.swap(other.m_bv);
        }

        void save(std::ostream& os) const
        {
            os.write(reinterpret_cast<char const*>(&m_n), sizeof(m_n));
            os.write(reinterpret_cast<char const*>(&m_hash_domain),
                     sizeof(m_hash_domain));
            m_hasher.save(os);
            m_bv.save(os);
        }

        void load(std::istream& is)
        {
            is.read(reinterpret_cast<char*>(&m_n), sizeof(m_n));
            is.read(reinterpret_cast<char*>(&m_hash_domain),
                    sizeof(m_hash_domain));
            m_hasher.load(is);
            m_bv.load(is);
        }


    private:

        uint64_t m_n;
        uint64_t m_hash_domain;
        BaseHasher m_hasher;
        ranked_bitpair_vector m_bv;
    };
}
