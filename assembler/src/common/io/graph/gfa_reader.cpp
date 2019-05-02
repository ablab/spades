//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gfa_reader.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/construction_helper.hpp"

#include "io/id_mapper.hpp"

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
                         io::IdMapper<std::string> *id_mapper) {
    auto helper = g.GetConstructionHelper();

    // INFO("Loading segments");
    std::vector<EdgeId> edges;
    edges.reserve(gfa_->n_seg);
    g.ereserve(2 * gfa_->n_seg);
    for (size_t i = 0; i < gfa_->n_seg; ++i) {
        gfa_seg_t *seg = gfa_->seg + i;

        uint8_t *kc = gfa_aux_get(seg->aux.l_aux, seg->aux.aux, "KC");
        unsigned cov = 0;
        if (kc && kc[0] == 'i')
            cov = *(int32_t*)(kc+1);
        DeBruijnEdgeData edata(Sequence(seg->seq));
        EdgeId e = helper.AddEdge(edata);
        g.coverage_index().SetRawCoverage(e, cov);
        g.coverage_index().SetRawCoverage(g.conjugate(e), cov);

        if (id_mapper) {
            (*id_mapper)[e.int_id()] = seg->name;
            if (e != g.conjugate(e)) {
                (*id_mapper)[g.conjugate(e).int_id()] = std::string(seg->name) + '\'';
            }
        }
        edges.push_back(e);
    }

    // INFO("Creating vertices");
    g.vreserve(gfa_->n_seg * 4);
    std::unordered_set<VertexId> vertices;
    for (uint32_t i = 0; i < gfa_->n_seg; ++i) {
        VertexId v1 = helper.CreateVertex(DeBruijnVertexData()),
                 v2 = helper.CreateVertex(DeBruijnVertexData());

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

    // INFO("Reading paths")
    paths_.reserve(gfa_->n_path);
    for (uint32_t i = 0; i < gfa_->n_path; ++i) {
        const gfa_path_t &path = gfa_->path[i];
        paths_.emplace_back(path.name);
        GFAPath &cpath = paths_.back();
        for (unsigned j = 0; j < path.n_seg; ++j) {
            EdgeId e = edges[path.v[j] >> 1];
            if (path.v[j] & 1)
                e = g.conjugate(e);
            cpath.edges.push_back(e);
        }
    }
}

}
