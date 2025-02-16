//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"

#include "assembly_graph/core/graph.hpp"
#include "sequence/sequence.hpp"

namespace io {

namespace binary {

template<typename Graph>
class GraphIO : public IOSingle<Graph> {
public:
    GraphIO()
            : IOSingle<Graph>("debruijn graph", ".grseq") {
    }

private:
    void SaveImpl(BinOStream &str, const Graph &graph) override {
        str << graph.vreserved() << graph.ereserved() << graph.link_size();

        size_t vertex_cnt = graph.size();
        str << vertex_cnt;

        std::vector<bool> saved_vertices(graph.max_vid(), false);

        auto SaveVertex = [&](typename Graph::VertexId v) {
          TRACE("Saving " << v.int_id() << " ~ " << graph.conjugate(v).int_id());
          str << v.int_id() << graph.conjugate(v).int_id();
          if (saved_vertices[v.int_id()])
              return;
          TRACE("Saving link info");
          bool complex = graph.is_complex(v);
          str << complex;
          if (!complex) {
              unsigned ovl = unsigned(graph.link_length(v, typename Graph::EdgeId(), typename Graph::EdgeId()));
              str << ovl;
          } else {
              str << graph.links(v).size();
              TRACE("Saving " << graph.links(v).size() << " links from " << v.int_id() << ", " << graph.conjugate(v).int_id());
              for (const auto &link_id: graph.links(v)) {
                  const auto &link = graph.link(link_id);
                  auto first_e = link.link.first;
                  auto second_e = link.link.second;
                  str << first_e << second_e << graph.conjugate(first_e) << graph.conjugate(second_e) << link.overlap;
              }
          }
          saved_vertices[v.int_id()] = true;
          saved_vertices[graph.conjugate(v).int_id()] = true;
        };

        for (auto v1 : graph) {
            SaveVertex(v1);
            for (auto e1 : graph.OutgoingEdges(v1)) {
                auto e2 = graph.conjugate(e1);
                if (e2 < e1)
                    continue;
                str << e1.int_id() << e2.int_id();
                //<< graph.EdgeEnd(e1).int_id() << graph.EdgeStart(e2).int_id() << graph.EdgeNucls(e1);

                SaveVertex(graph.EdgeEnd(e1));
                str << graph.EdgeNucls(e1);
            }
            str << (size_t)0; //null-term
        }
    }

    void LoadImpl(BinIStream &str, Graph &graph) override {
        graph.clear();

        uint64_t max_vid, max_eid, num_links;
        str >> max_vid >> max_eid >> num_links;
        graph.reserve(max_vid, max_eid);
        TRACE("Reserving " << num_links << " links");
        graph.lreserve(num_links);

        size_t vertex_cnt;
        str >> vertex_cnt;

        auto TryAddVertex = [&](uint64_t ids[2]) {
            if (graph.contains(typename Graph::VertexId(ids[0])))
                return;
            TRACE("Vertex " << ids[0] << " ~ " << ids[1] << " .");
            bool complex;
            typename Graph::VertexId new_id;
            str >> complex;
            TRACE("Complex: " << complex);
            if (!complex) {
                unsigned ovl;
                str >> ovl;
                new_id = graph.AddVertex(typename Graph::VertexData(ovl), ids[0], ids[1]);
            } else {
                uint link_count = 0;
                str >> link_count;
                std::vector<debruijn_graph::LinkId> link_ids;
                std::vector<debruijn_graph::LinkId> conj_link_ids;
                TRACE("Reading " << link_count << " links from " << ids[0] << ", " << ids[1]);
                for (uint i = 0; i < link_count; ++i) {
                    typename Graph::EdgeId e1, e2, e1_conj, e2_conj;
                    unsigned ovl;
                    str >> e1 >> e2 >> e1_conj >> e2_conj >> ovl;
                    auto link_id = graph.add_link(e1, e2, ovl);
                    auto conj_link_id = graph.add_link(e2_conj, e1_conj, ovl);
                    link_ids.push_back(link_id);
                    conj_link_ids.push_back(conj_link_id);
                }
                std::vector<debruijn_graph::LinkId> empty_links;
                new_id = graph.AddVertex(debruijn_graph::DeBruijnVertexData(empty_links), ids[0], ids[1]);
                auto conj_id = graph.conjugate(new_id);
                graph.add_links(new_id, link_ids);
                graph.add_links(conj_id, conj_link_ids);
            }
            VERIFY(new_id == ids[0]);
            VERIFY(graph.conjugate(new_id) == ids[1]);
            TRACE("Added " << ids[0] << " ~ " << ids[1] << " .")
        };

        for (size_t i = 0; i < vertex_cnt; ++i) {
            uint64_t start_ids[2];
            str >> start_ids;
            TryAddVertex(start_ids);
            while (true) {
                uint64_t edge_ids[2];
                str >> edge_ids[0];
                if (!edge_ids[0]) //null-term
                    break;
                str >> edge_ids[1];
                uint64_t end_ids[2];
                Sequence seq;
                str >> end_ids;
                TRACE("Edge " << edge_ids[0] << " : " << start_ids[0] << " -> "
                              << end_ids[0] << " l = " << seq.size() << " ~ " << edge_ids[1]);
                TryAddVertex(end_ids);
                str >> seq;

                auto new_id = graph.AddEdge(start_ids[0], end_ids[0],
                        typename Graph::EdgeData(seq), edge_ids[0], edge_ids[1]);
                VERIFY(new_id == edge_ids[0]);
                VERIFY(graph.conjugate(new_id) == edge_ids[1]);
            }
        }
    }

    DECL_LOGGER("GraphIO");
};

template<>
struct IOTraits<debruijn_graph::Graph> {
    typedef GraphIO<debruijn_graph::Graph> Type;
};

} // namespace binary

} // namespace io
