#include "multiplex_gfa_reader.hpp"

#include "assembly_graph/core/construction_helper.hpp"

#include <unordered_set>

using namespace debruijn_graph;

namespace cont_index {

MultiplexGFAReader::MultiplexGFAReader()
    : gfa_(nullptr, gfa_destroy) {}
MultiplexGFAReader::MultiplexGFAReader(const std::string &filename)
    : gfa_(gfa_read(filename.c_str()), gfa_destroy) {}
bool MultiplexGFAReader::open(const std::string &filename) {
    gfa_.reset(gfa_read(filename.c_str()));

    return (bool)gfa_;
}

uint32_t MultiplexGFAReader::num_edges() const { return gfa_->n_seg; }
uint64_t MultiplexGFAReader::num_links() const { return gfa_->n_arc; }

unsigned MultiplexGFAReader::k() const {
    unsigned k = -1U;
    for (size_t i = 0; i < gfa_->n_arc; ++i) {
        gfa_arc_t *arc = gfa_->arc + i;
        if (arc->ov != arc->ow || arc->ov < 0)
            return -1U;

        if (k == -1U || k > unsigned(arc->ov))
            k = unsigned(arc->ov);
    }

    return k;
}

void MultiplexGFAReader::to_graph(ConjugateDeBruijnGraph &g,
                                  io::IdMapper<std::string> *id_mapper) {
    auto helper = g.GetConstructionHelper();

     INFO("Loading segments");
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
            DEBUG("Map ids: " << e.int_id() << ":" << seg->name);
            DEBUG("Map ids: " << g.conjugate(e).int_id() << ":" << seg->name << "'");
            (*id_mapper)[e.int_id()] = seg->name;
            if (e != g.conjugate(e)) {
                (*id_mapper)[g.conjugate(e).int_id()] = std::string(seg->name) + '\'';
            }
        }
        edges.push_back(e);
    }

     INFO("Creating vertices");
    g.vreserve(gfa_->n_seg * 4);
    std::unordered_set<VertexId> vertices;
    for (uint32_t i = 0; i < gfa_->n_seg; ++i) {
        VertexId v1 = helper.CreateVertex(DeBruijnVertexData());
        helper.LinkIncomingEdge(v1, edges[i]);
        vertices.insert(v1);

        if (edges[i] != g.conjugate(edges[i])) {
            VertexId v2 = helper.CreateVertex(DeBruijnVertexData());
            helper.LinkIncomingEdge(v2, g.conjugate(edges[i]));
            vertices.insert(v2);
        }
    }

    INFO("Linking edges");
    std::unordered_map<VertexId, size_t> vertex_to_overlap;
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

                VERIFY(av[j].ov == av[j].ow && av[j].ov >= g.k());
                size_t overlap = static_cast<size_t>(av[j].ov);
                VertexId link_vertex = g.EdgeEnd(e1);
                auto overlap_result = vertex_to_overlap.find(link_vertex);
                if (overlap_result == vertex_to_overlap.end()) {
                    vertex_to_overlap.emplace(link_vertex, overlap);
                } else {
                    VERIFY(overlap_result->second == overlap);
                }
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

                VERIFY(av[j].ov == av[j].ow && av[j].ov >= g.k());
                size_t overlap = static_cast<size_t>(av[j].ov);
                VertexId link_vertex = g.EdgeEnd(e1);
                auto overlap_result = vertex_to_overlap.find(link_vertex);
                if (overlap_result == vertex_to_overlap.end()) {
                    vertex_to_overlap.emplace(link_vertex, overlap);
                } else {
                    VERIFY(overlap_result->second == overlap);
                }
            }
        }
    }

    //fixme just add this to regular GFAReader
    for (const auto &vertex: vertices) {
        if (g.OutgoingEdgeCount(vertex) > 0 and g.IncomingEdgeCount(vertex) > 0) {
            size_t overlap = vertex_to_overlap.at(vertex);
            auto first_incoming = *(g.IncomingEdges(vertex).begin());
            Sequence overlap_sequence = g.EdgeNucls(first_incoming).Last(overlap);
            DeBruijnEdgeData edata(overlap_sequence);
            EdgeId overlap_edge = helper.AddEdge(edata);
            double incoming_cov = .0;
            double outgoing_cov = .0;

            if (g.conjugate(overlap_edge) == overlap_edge) {
                VERIFY_MSG(false, "Self-conjugate edges are currently not supported");
            }
            VertexId overlap_end = helper.CreateVertex(DeBruijnVertexData());
            helper.LinkOutgoingEdge(vertex, overlap_edge);
            helper.LinkIncomingEdge(overlap_end, overlap_edge);

            //fixme ineffective
            for (const auto &incoming: g.IncomingEdges(vertex)) {
                Sequence current_overlap = g.EdgeNucls(incoming).Last(overlap);
                VERIFY_MSG(current_overlap == overlap_sequence, "Input graph is not multiplex DBG!");
                incoming_cov += g.coverage(incoming);
            }
            for (const auto &outgoing: g.OutgoingEdges(vertex)) {
                Sequence current_overlap = g.EdgeNucls(outgoing).First(overlap);
                VERIFY_MSG(current_overlap == overlap_sequence, "Input graph is not multiplex DBG!");
                outgoing_cov += g.coverage(outgoing);
                helper.LinkEdges(overlap_edge, outgoing);
            }
            double cov = (incoming_cov + outgoing_cov) / 2;

            g.coverage_index().SetRawCoverage(overlap_edge, cov);
            g.coverage_index().SetRawCoverage(g.conjugate(overlap_edge), cov);
        }
    }

    // INFO("Filtering dangling vertices");
    for (VertexId v : vertices) {
        if (g.OutgoingEdgeCount(v) > 0 || g.IncomingEdgeCount(v) > 0)
            continue;

        g.DeleteVertex(v);
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
