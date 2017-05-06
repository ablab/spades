//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/stats/picture_dump.hpp"
#include <io/reads/osequencestream.hpp>
#include "assembly_graph/components/connected_component.hpp"
#include "assembly_graph/stats/statistics.hpp"
#include "assembly_graph/paths/path_finders.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/paths/bidirectional_path_io/io_support.hpp"

namespace debruijn_graph {

template<class Graph>
class GFAWriter {
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    std::ostream &os_;

    bool IsCanonical(EdgeId e) const {
        return e <= graph_.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : graph_.conjugate(e);
    }

    std::string GetOrientation(EdgeId e) const {
        return IsCanonical(e) ? "+" : "-";
    }

    void WriteSegment(size_t edge_id, const Sequence &seq, double cov) {
        os_ << "S\t" << edge_id << "\t"
           << seq.str() << "\t"
           << "KC:i:" << size_t(math::round(cov)) << std::endl;
    }

    void WriteSegments() {
        for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            WriteSegment(graph_.int_id(e), graph_.EdgeNucls(e), graph_.coverage(e) * double(graph_.length(e)));
        }
    }

    void WriteLink(EdgeId e1, EdgeId e2,
                   size_t overlap_size) {
        os_ << "L\t" << EdgeString(e1, "\t") << "\t"
           << EdgeString(e2, "\t") << "\t"
           << overlap_size << "M" << std::endl;
    }

    void WriteLinks() {
        for (VertexId v : graph_)
            for (auto inc_edge : graph_.IncomingEdges(v))
                for (auto out_edge : graph_.OutgoingEdges(v))
                    WriteLink(inc_edge, out_edge, graph_.k());
    }

    std::string EdgeString(EdgeId e, const std::string &delim = "") {
        return std::to_string(Canonical(e).int_id()) + delim + GetOrientation(e);
    }

    void WritePath(size_t path_id, size_t segment_id, const vector<std::string> &edge_strs) {
        os_ << "P" << "\t" ;
        os_ << path_id << "_" << segment_id << "\t";
        std::string delimeter = "";
        for (const auto& e : edge_strs) {
            os_ << delimeter << e;
            delimeter = ",";
        }
        os_ << "\t";
        delimeter = "";
        for (size_t i = 0; i < edge_strs.size() - 1; ++i) {
            os_ << delimeter << "*";
            delimeter = ",";
        }
        os_ << std::endl;
    }

public:
    GFAWriter(const Graph &graph, std::ostream &os)
            : graph_(graph), os_(os) {
    }

    void WriteSegmentsAndLinks() {
        WriteSegments();
        WriteLinks();
    }

    void WritePaths(const path_extend::PathContainer &paths) {
        for (const auto &path_pair : paths) {
            const path_extend::BidirectionalPath &p = (*path_pair.first);
            if (p.Size() == 0) {
                continue;
            }
            std::vector<std::string> segmented_path;
            size_t id = p.GetId();
            size_t segment_id = 1;
            for (size_t i = 0; i < p.Size() - 1; ++i) {
                EdgeId e = p[i];
                segmented_path.push_back(EdgeString(e));
                if (graph_.EdgeEnd(e) != graph_.EdgeStart(p[i+1]) || p.GapAt(i+1).gap > 0) {
                    WritePath(id, segment_id, segmented_path);
                    segment_id++;
                    segmented_path.clear();
                }
            }

            segmented_path.push_back(EdgeString(p.Back()));
            WritePath(id, segment_id, segmented_path);
        }
    }

};

inline void OutputContigs(const Graph &g,
                          const string &contigs_output_filename) {
    INFO("Outputting contigs to " << contigs_output_filename << ".fasta");
    io::osequencestream_cov oss(contigs_output_filename + ".fasta");

    for (auto it = g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        oss << g.EdgeNucls(e).str();
        oss << g.coverage(e);
    }
}

}
