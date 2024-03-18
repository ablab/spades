//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "link_index.hpp"
#include "io/utils/id_mapper.hpp"

using namespace binning;

void GraphLinkIndex::Init(const debruijn_graph::Graph &g) {
    for (EdgeId e : g.canonical_edges()) {
        for (EdgeId o : g_.OutgoingEdges(g.EdgeEnd(e)))
            add(e, o);

        for (EdgeId i : g_.IncomingEdges(g.EdgeStart(e)))
            add(e, i);
    }
}

void LinkIndex::dump(const std::string &output_path, const io::IdMapper<std::string> &edge_mapper) {
    std::ofstream os(output_path);

    os << "FirstId\tSecondId\tWeight\tFirstLength\tFirstCov\tSecondLength\tSecondCov\n";
    for (auto it = data_.cbegin(), end = data_.cend(); it != end; ++it) {
        EdgeId edge = it.key();

        for (const auto &neighbour: it.value()) {
            EdgeId neigh = neighbour.e;
            os << edge_mapper[edge.int_id()] << "\t" << edge_mapper[neigh.int_id()] << "\t" << neighbour.w << "\t"
               << g_.length(edge) + g_.k() - 1 << "\t" << g_.coverage(edge) << "\t"
               << g_.length(neigh) + g_.k() - 1 << "\t" << g_.coverage(neigh)  << "\n";
        }
    }
}
