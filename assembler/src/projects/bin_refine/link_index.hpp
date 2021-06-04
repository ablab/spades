//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "id_map.hpp"
#include "assembly_graph/core/graph.hpp"

#include "adt/small_pod_vector.hpp"
#include "adt/flat_set.hpp"

namespace io {
template<class T>
class IdMapper;
}

namespace binning {

class LinkIndex {
  public:
    using EdgeId = debruijn_graph::EdgeId;

    struct EdgeWithWeight {
        EdgeId e;
        double w;

        EdgeWithWeight() = default;
        EdgeWithWeight(EdgeId e_, double w_)
                : e(e_), w(w_) {}

        bool operator<(const EdgeWithWeight &rhs) const {
            return e < rhs.e;
        }
    };
    using EdgeLinks = adt::flat_set<EdgeWithWeight, std::less<EdgeWithWeight>, adt::SmallPODVector>;

    LinkIndex(const debruijn_graph::Graph &g)
            : g_(g), data_(g.max_eid()) {}

    void add(EdgeId e1, EdgeId e2, double w = 1.) {
        // Link index must be symmetric and all links must be unique
        data_[e1].emplace(EdgeWithWeight{e2, w});
        data_[e2].emplace(EdgeWithWeight{e1, w});
    }

    const EdgeLinks &links(EdgeId e1) const {
        return data_.at(e1);
    }

    void dump(const std::string &output_path, const io::IdMapper<std::string> &edge_mapper);

  protected:
    const debruijn_graph::Graph &g_;
  private:
    adt::id_map<EdgeLinks, EdgeId> data_;
};


class GraphLinkIndex : public LinkIndex {
  public:
    using LinkIndex::EdgeId;

    GraphLinkIndex(const debruijn_graph::Graph &g)
            : LinkIndex(g) {
        Init(g);
    }

    void Init(const debruijn_graph::Graph &g);

    void add(EdgeId e1, EdgeId e2, double w = 1.) {
        LinkIndex::add(e1, e2, w);
        EdgeId ce1 = g_.conjugate(e1);
        if (e1 != ce1)
            LinkIndex::add(ce1, g_.conjugate(e2));
    }
};


}
