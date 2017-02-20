#pragma once

#include <cassert>
#include <cstdint>
#include <tuple>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include "common.hpp"
#include "hypergraph.hpp"

#include "utils/logger/logger.hpp"

namespace emphf {

    template <typename HypergraphType>
    class hypergraph_sorter_seq {
    public:
        typedef HypergraphType hg;
        typedef typename hg::node_t node_t;
        typedef typename hg::hyperedge hyperedge;
        typedef typename hg::xored_adj_list xored_adj_list;

        hypergraph_sorter_seq()
        {}

        template <typename Range, typename EdgeGenerator>
        bool try_generate_and_sort(Range const& input_range,
                                   EdgeGenerator const& edge_gen,
                                   size_t n,
                                   size_t hash_domain,
                                   bool verbose = true)
        {
            using std::get;
            std::vector<xored_adj_list> adj_lists;

            size_t m = hash_domain * 3;

            // do all the allocations upfront
            m_peeling_order.clear();
            m_peeling_order.reserve(n);
            adj_lists.resize(m);

            // generate edges
            if (verbose) {
                //logger() << "Generating hyperedges and populating adjacency lists"
                //         << std::endl;
            }

            for (auto const& val: input_range) {
                auto edge = edge_gen(val);
                // canonical by construction
                assert(orientation(edge) == 0);

                adj_lists[edge.v0].add_edge(edge);

                std::swap(edge.v0, edge.v1);
                adj_lists[edge.v0].add_edge(edge);

                std::swap(edge.v0, edge.v2);
                adj_lists[edge.v0].add_edge(edge);
            }

            // peel
            if (verbose) {
                // logger() << "Peeling" << std::endl;
            }

            auto visit = [&](node_t v0) {
                if (adj_lists[v0].degree == 1) {
                    auto edge = adj_lists[v0].edge_from(v0);
                    m_peeling_order.push_back(edge);

                    edge = canonicalize_edge(edge);
                    adj_lists[edge.v0].delete_edge(edge);

                    std::swap(edge.v0, edge.v1);
                    adj_lists[edge.v0].delete_edge(edge);

                    std::swap(edge.v0, edge.v2);
                    adj_lists[edge.v0].delete_edge(edge);
                }
            };

            size_t queue_position = 0;
            for (node_t v0 = 0; v0 < m; ++v0) {
                visit(v0);

                while (queue_position < m_peeling_order.size()) {
                    auto const& cur_edge = m_peeling_order[queue_position];

                    visit(cur_edge.v1);
                    visit(cur_edge.v2);
                    queue_position += 1;
                }
            }

            if (m_peeling_order.size() < n) {
                if (verbose) {
                    // logger() << "Hypergraph is not peelable: "
                    //         << (n - m_peeling_order.size()) << " edges remaining"
                    //         << std::endl;
                }
                return false;
            }

            assert(m_peeling_order.size() == n);

            return true;
        }

        typedef typename std::vector<hyperedge>::const_reverse_iterator
        peeling_iterator;

        std::pair<peeling_iterator, peeling_iterator>
        get_peeling_order() const
        {
            return std::make_pair(m_peeling_order.crbegin(),
                                  m_peeling_order.crend());
        }

    private:

        size_t m_hash_domain;
        std::vector<hyperedge> m_peeling_order;
    };
}
