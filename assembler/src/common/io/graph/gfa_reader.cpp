//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gfa_reader.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/construction_helper.hpp"

#include "gfa1/gfa.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_set>

using namespace debruijn_graph;

namespace gfa {

GFAReader::GFAReader()
        : gfa_(nullptr, gfa_destroy) {}
GFAReader::GFAReader(const std::string &filename)
        : gfa_(gfa_read(filename.c_str()), gfa_destroy) {}
bool GFAReader::open(const std::string &filename) {
    gfa_.reset(gfa_read(filename.c_str()));

    return (bool)gfa_;
}

uint32_t GFAReader::num_edges() const { return gfa_->n_seg; }
uint64_t GFAReader::num_links() const { return gfa_->n_arc; }

void GFAReader::to_graph(ConjugateDeBruijnGraph &g,
                         bool numeric_ids) {
    auto helper = g.GetConstructionHelper();

    // INFO("Loading segments");
    std::vector<EdgeId> edges;
    edges.reserve(gfa_->n_seg);

    if (numeric_ids) {
        uint64_t mid = 0;
        std::unordered_set<uint64_t> ids;
        for (size_t i = 0; i < gfa_->n_seg; ++i) {
            gfa_seg_t *seg = gfa_->seg + i;
            uint64_t id = std::atoll(seg->name);
            VERIFY_MSG(ids.insert(id).second, "Unique numeric ids are required");
            if (mid < id)
                mid = id;
        }
        restricted::IdSegmentStorage eid_storage = helper.graph().GetGraphIdDistributor().Reserve(mid + 2, 
                /*force zero shift*/true);
        for (size_t i = 0; i < gfa_->n_seg; ++i) {
            gfa_seg_t *seg = gfa_->seg + i;

            uint64_t id = std::atoll(seg->name);

            uint64_t ids[] = { id, id + 1};
            auto id_distributor = eid_storage.GetSegmentIdDistributor(std::begin(ids), std::end(ids));
            EdgeId e = helper.AddEdge(DeBruijnEdgeData(Sequence(seg->seq)), id_distributor);
            edges.push_back(e);
        }
    } else {
        restricted::IdSegmentStorage eid_storage = helper.graph().GetGraphIdDistributor().Reserve(gfa_->n_seg * 2);
        for (size_t i = 0; i < gfa_->n_seg; ++i) {
            gfa_seg_t *seg = gfa_->seg + i;

            auto id_distributor = eid_storage.GetSegmentIdDistributor(i << 1, (i << 1) + 2);
            EdgeId e = helper.AddEdge(DeBruijnEdgeData(Sequence(seg->seq)), id_distributor);
            edges.push_back(e);
        }
    }

    // INFO("Creating vertices");
    restricted::IdSegmentStorage vid_storage = helper.graph().GetGraphIdDistributor().Reserve(gfa_->n_seg * 4);
    std::unordered_set<VertexId> vertices;
    for (uint32_t i = 0; i < gfa_->n_seg; ++i) {
        auto id_distributor = vid_storage.GetSegmentIdDistributor(size_t(i << 2), size_t(i << 2) + 4); // indices for four vertices are required
        VertexId v1 = helper.CreateVertex(DeBruijnVertexData(), id_distributor),
                 v2 = helper.CreateVertex(DeBruijnVertexData(), id_distributor);

        helper.LinkIncomingEdge(v1, edges[i]);
        if (edges[i] != g.conjugate(edges[i]))
            helper.LinkIncomingEdge(v2, g.conjugate(edges[i]));

        vertices.insert(v1);
        vertices.insert(v2);
    }

    // INFO("Linking edges");
    for (uint32_t i = 0; i < gfa_->n_seg; ++i) {
        EdgeId e1 = edges[i];
        // Process direct links
        {
            uint32_t vv = i << 1 | 0;
            gfa_arc_t *av = gfa_arc_a(gfa_.get(), vv);
            for (size_t j = 0; j < gfa_arc_n(gfa_.get(), vv); ++j) {
                EdgeId e2 = edges[av[j].w >> 1];
                if (av[j].w & 1)
                    e2 = g.conjugate(e2);
                helper.LinkEdges(e1, e2);
            }
        }

        // Process rc links
        {
            e1 = g.conjugate(e1);
            uint32_t vv = i << 1 | 1;
            gfa_arc_t *av = gfa_arc_a(gfa_.get(), vv);
            for (size_t j = 0; j < gfa_arc_n(gfa_.get(), vv); ++j) {
                EdgeId e2 = edges[av[j].w >> 1];
                if (av[j].w & 1)
                    e2 = g.conjugate(e2);
                helper.LinkEdges(e1, e2);
            }
        }
    }

    // INFO("Filtering dangling vertices");
    for (auto it = vertices.begin(); it != vertices.end(); ) {
        VertexId v = *it;
        if (g.OutgoingEdgeCount(v) == 0 && g.IncomingEdgeCount(v) == 0) {
            it = vertices.erase(it);
        } else
            ++it;
    }
    helper.AddVerticesToGraph(vertices.begin(), vertices.end());
}

}
