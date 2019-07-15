#include "contracted_gfa_writer.hpp"

namespace path_extend {
namespace read_cloud {
void ContractedGFAWriter::WriteLinks() {
    for (const auto &vertex: graph_) {
        for (const auto &in_entry: graph_.incoming(vertex)) {
            for (const auto &in_edge: in_entry.second) {
                if (edge_namer_.IsCanonical(in_edge)) {
                    INFO("Writing link");
                    WriteLink(in_edge, vertex, graph_.GetAssemblyGraph().k(), os_, edge_namer_);
                }
            }
        }
        for (const auto &out_entry: graph_.outcoming(vertex)) {
            for (const auto &out_edge: out_entry.second) {
                if (edge_namer_.IsCanonical(out_edge)) {
                    WriteLink(vertex, out_edge, graph_.GetAssemblyGraph().k(), os_, edge_namer_);
                }
            }
        }
    }
}
void ContractedGFAWriter::WriteLink(ContractedGFAWriter::EdgeId edge,
                                    VertexId vertex,
                                    size_t overlap_size,
                                    std::ostream &os,
                                    io::CanonicalEdgeHelper<ContractedGFAWriter::ContractedGraph> &namer) {
    string vertex_orientation_string = "NODE_" + std::to_string(vertex.int_id()) + "\t+";
    os << "L\t"
       << namer.EdgeOrientationString(edge, "\t") << '\t'
       << vertex_orientation_string << '\t'
       << overlap_size << "M\n";
}
void ContractedGFAWriter::WriteLink(VertexId vertex,
                                    ContractedGFAWriter::EdgeId edge,
                                    size_t overlap_size,
                                    std::ostream &os,
                                    io::CanonicalEdgeHelper<ContractedGFAWriter::ContractedGraph> &namer) {
    string vertex_orientation_string = "NODE_" + std::to_string(vertex.int_id()) + "\t+";
    os << "L\t"
       << vertex_orientation_string << '\t'
       << namer.EdgeOrientationString(edge, "\t") << '\t'
       << overlap_size << "M\n";
}
void ContractedGFAWriter::WriteSegments() {
    for (const auto &vertex: graph_) {
        if (graph_.GetInDegree(vertex) > 0 or graph_.GetOutDegree(vertex) > 0) {
            std::string vertex_sequence = graph_.GetAssemblyGraph().VertexNucls(vertex).str();
            const double vertex_coverage = 100.0;
            string vertex_name = "NODE_" + std::to_string(vertex.int_id());
            WriteSegment(vertex_name, vertex_sequence, vertex_coverage, os_);
        }
    }
    for (const auto &vertex: graph_) {
        for (const auto &out_entry: graph_.outcoming(vertex)) {
            for (const auto &edge: out_entry.second) {
                std::string sequence = graph_.EdgeNucls(edge);
                if (edge_namer_.IsCanonical(edge) and sequence.length() > 0) {
                    WriteSegment(edge_namer_.EdgeString(edge), graph_.EdgeNucls(edge),
                                 graph_.coverage(edge) * double(graph_.length(edge)),
                                 os_);
                }
            }
        }
    }
}
void ContractedGFAWriter::WriteSegment(const std::string &edge_id,
                                       const std::string &seq,
                                       double cov,
                                       std::ostream &os) {
    os << "S\t"
       << edge_id << '\t' << seq << '\t'
       << "KC:i:" << size_t(math::round(cov)) << '\n';
}
void ContractedGFAWriter::WriteSegmentsAndLinks() {
    WriteSegments();
    WriteLinks();
}
//void ContractedCSVWriter::WriteProperties() {
//    for (const auto &vertex: graph_) {
//        graph_.GetAssemblyGraph().conjugate(vertex);
//    }
//}
}
}