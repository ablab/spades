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

void MultiplexGFAReader::SafeInsert(OverlapStorage &storage, const std::string &name, const Sequence &seq) {
    auto result = storage.find(name);
    if (result == storage.end()) {
        storage.insert({name, seq});
    } else {
        VERIFY_MSG(result->second == seq, "Inconsistency in segment " << name << ", " << seq.size() << ", " << result->second.size());
    }
}

MultiplexGFAReader::OverlapStorages MultiplexGFAReader::ConstructOverlapStorages() const {
    OverlapStorage start_to_seq;
    OverlapStorage end_to_seq;

    for (uint32_t i = 0; i < gfa_->n_seg; ++i) {
        {
            uint32_t vv = i << 1 | 0;
            gfa_arc_t *av = gfa_arc_a(gfa_.get(), vv);
            gfa_seg_t *from = gfa_->seg + i;
            for (size_t j = 0; j < gfa_arc_n(gfa_.get(), vv); ++j) {
                size_t overlap = static_cast<size_t>(av[j].ov);

                std::string from_name = from->name;
//                INFO("From " << from_name);
                Sequence from_seq = Sequence(from->seq).Last(overlap);
                SafeInsert(end_to_seq, from_name, from_seq);
                gfa_seg_t *to = gfa_->seg + (av[j].w >> 1);
                std::string to_name = av[j].w & 1 ? std::string(to->name) + "\'" : to->name;
//                INFO("To " << to_name);
                Sequence to_seq = av[j].w & 1 ? !(Sequence(to->seq).Last(overlap)) : Sequence(to->seq).First(overlap);
                SafeInsert(start_to_seq, to_name, to_seq);
            }
        }

        {
            uint32_t vv = i << 1 | 1;
            gfa_arc_t *av = gfa_arc_a(gfa_.get(), vv);
            gfa_seg_t *from = gfa_->seg + i;
            for (size_t j = 0; j < gfa_arc_n(gfa_.get(), vv); ++j) {
                size_t overlap = static_cast<size_t>(av[j].ov);

                std::string from_name = std::string(from->name) + "\'";
//                INFO("From " << from_name);
                Sequence from_seq = !(Sequence(from->seq).First(overlap));
                SafeInsert(end_to_seq, from_name, from_seq);
                gfa_seg_t *to = gfa_->seg + (av[j].w >> 1);
                std::string to_name = av[j].w & 1 ? std::string(to->name) + "\'" : to->name;
//                INFO("To " << to_name);
                Sequence to_seq = av[j].w & 1 ? !(Sequence(to->seq).Last(overlap)) : Sequence(to->seq).First(overlap);
                SafeInsert(start_to_seq, to_name, to_seq);
            }
        }
    }
    OverlapStorages result(start_to_seq, end_to_seq);
    return result;
}

void MultiplexGFAReader::to_graph(ConjugateDeBruijnGraph &g,
                                  io::IdMapper<std::string> *id_mapper) {
    auto helper = g.GetConstructionHelper();

    auto overlap_storages = ConstructOverlapStorages();
    const auto &start_overlap_storage = overlap_storages.first;
    const auto &end_overlap_storage = overlap_storages.second;
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

        size_t start_overlap = 0;
        size_t end_overlap = 0;
        auto start_ov_result = start_overlap_storage.find(seg->name);
        auto end_ov_result = end_overlap_storage.find(seg->name);
        if (start_ov_result != start_overlap_storage.end()) {
            start_overlap = start_ov_result->second.size();
        }
        if (end_ov_result != end_overlap_storage.end()) {
            end_overlap = end_ov_result->second.size();
        }
        if (start_overlap > 0) {
            VERIFY(start_overlap >= g.k());
        } else {
            start_overlap = g.k();
        }
        if (end_overlap > 0) {
            VERIFY(end_overlap >= g.k());
        } else {
            end_overlap = g.k();
        }
        INFO(seg->name);
        INFO("Overlaps: " << start_overlap << ", " << end_overlap);
        auto init_seq = Sequence(seg->seq);
        size_t init_size = init_seq.size();
        size_t from = start_overlap - g.k();
        size_t to = init_size - end_overlap + g.k() - 1;
        INFO(init_size << ", " << from << ", " << to);
        DeBruijnEdgeData edata(Sequence(seg->seq).Subseq(from, to));

        EdgeId e = helper.AddEdge(edata);
        g.coverage_index().SetRawCoverage(e, cov);
        g.coverage_index().SetRawCoverage(g.conjugate(e), cov);
        INFO("Constructed edge");

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

    std::vector<VertexId> vertices_copy;
    std::copy(g.vertices().begin(), g.vertices().end(), std::back_inserter(vertices_copy));
    for (const auto &vertex: vertices_copy) {
        if (g.OutgoingEdgeCount(vertex) > 0 and g.IncomingEdgeCount(vertex) > 0) {
            auto first_incoming = *(g.IncomingEdges(vertex).begin());
            bool is_canonical = first_incoming <= g.conjugate(first_incoming);
            auto first_name = is_canonical ? (*id_mapper)[first_incoming.int_id()] : (*id_mapper)[first_incoming.int_id()] + "\'";
            Sequence overlap_sequence = end_overlap_storage.at(first_name);
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

            double cov = (incoming_cov + outgoing_cov) / 2;

            g.coverage_index().SetAvgCoverage(overlap_edge, cov);
            g.coverage_index().SetAvgCoverage(g.conjugate(overlap_edge), cov);
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
