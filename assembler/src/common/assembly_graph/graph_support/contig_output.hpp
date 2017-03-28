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

namespace debruijn_graph {

class GFASegmentWriter {
private:
    std::ostream &ostream_;

public:

    GFASegmentWriter(std::ostream &stream) : ostream_(stream)  {
    }

    void Write(size_t edge_id, const Sequence &seq, double cov) {
        ostream_ << "S\t" << edge_id << "\t";
        ostream_ << seq.str() << "\t";
        ostream_ << "KC:i:" << int(cov) << std::endl;
    }
};

class GFALinkWriter {
private:
    std::ostream &ostream_;
    size_t overlap_size_;

public:

    GFALinkWriter(std::ostream &stream, size_t overlap_size) : ostream_(stream), overlap_size_(overlap_size)  {
    }

    void Write(size_t first_segment, std::string &first_orientation, size_t second_segment, std::string &second_orientation) {
        ostream_ << "L\t" << first_segment << "\t" << first_orientation << "\t" ;
        ostream_ << second_segment << "\t" << second_orientation << "\t" << overlap_size_ << "M";
        ostream_ << std::endl;

    }
};

struct PathSegmentSequence {
    size_t path_id_;
    size_t segment_number_;
    std::vector<std::string> segment_sequence_;
    PathSegmentSequence(size_t path_id, std::vector<std::string> &segment_sequence)
    : path_id_(path_id), segment_number_(1), segment_sequence_(segment_sequence) {
    }

    PathSegmentSequence()
    : path_id_(0), segment_number_(1), segment_sequence_(){
    }
    void Reset() {
        segment_sequence_.clear();
    }
};

class GFAPathWriter {
    std::ostream &ostream_;

public:

    GFAPathWriter(std::ostream &stream)
    : ostream_(stream)  {
    }

    void Write(const PathSegmentSequence &path_segment_sequence) {
        ostream_ << "P" << "\t" ;
        ostream_ << path_segment_sequence.path_id_ << "_" << path_segment_sequence.segment_number_ << "\t";
        std::string delimeter = "";
        for (size_t i = 0; i < path_segment_sequence.segment_sequence_.size(); ++i) {
            ostream_ << delimeter << path_segment_sequence.segment_sequence_[i];
            delimeter = ",";
        }
        ostream_ << "\t";
        delimeter = "";
        for (size_t i = 0; i < path_segment_sequence.segment_sequence_.size() - 1; ++i) {
            ostream_ << delimeter << "*";
            delimeter = ",";
        }
        ostream_ << std::endl;
    }
};

template<class Graph>
class GFAWriter {
private:
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    const path_extend::PathContainer &paths_;
    const string filename_;

    bool IsCanonical(EdgeId e) const {
        if (e <= graph_.conjugate(e)) {
            return true;
        } else {
            return false;
        }
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : graph_.conjugate(e);
    }

    std::string GetOrientation(EdgeId e) const {
        return IsCanonical(e) ? "+" : "-";
    }

    void WriteSegments(std::ofstream &stream) {
        GFASegmentWriter segment_writer(stream);
        for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            segment_writer.Write((*it).int_id(), graph_.EdgeNucls(*it), graph_.coverage(*it) * double(graph_.length(*it)));
        }
    }

    void WriteLinks(std::ofstream &stream) {
        GFALinkWriter link_writer(stream, graph_.k());
        for (auto it = graph_.SmartVertexBegin(); !it.IsEnd(); ++it) {
            for (auto inc_edge : graph_.IncomingEdges(*it)) {
                std::string orientation_first = GetOrientation(inc_edge);
                size_t segment_first = Canonical(inc_edge).int_id();
                for (auto out_edge : graph_.OutgoingEdges(*it)) {
                    size_t segment_second = Canonical(out_edge).int_id();
                    std::string orientation_second = GetOrientation(out_edge);
                    link_writer.Write(segment_first, orientation_first, segment_second, orientation_second);
                }
            }
        }
    }

    void UpdateSegmentedPath(PathSegmentSequence &segmented_path, EdgeId e) {
        std::string segment_id = std::to_string(Canonical(e).int_id());
        std::string orientation = GetOrientation(e);
        segmented_path.segment_sequence_.push_back(segment_id + orientation);
    }

    void WritePaths(std::ofstream &stream) {
        GFAPathWriter path_writer(stream);
        for (const auto &path_pair : paths_) {
            const path_extend::BidirectionalPath &p = (*path_pair.first);
            if (p.Size() == 0) {
                continue;
            }
            PathSegmentSequence segmented_path;
            segmented_path.path_id_ = p.GetId();
            for (size_t i = 0; i < p.Size() - 1; ++i) {
                EdgeId e = p[i];
                UpdateSegmentedPath(segmented_path, e);
                if (graph_.EdgeEnd(e) != graph_.EdgeStart(p[i+1])) {
                    path_writer.Write(segmented_path);
                    segmented_path.segment_number_++;
                    segmented_path.Reset();
                }
            }
            UpdateSegmentedPath(segmented_path, p.Back());
            path_writer.Write(segmented_path);

        }
    }

public:
    GFAWriter(const Graph &graph, const path_extend::PathContainer &paths, const string &filename)
    : graph_(graph), paths_(paths), filename_(filename) {
    }

    void Write() {
        std::ofstream stream;
        stream.open(filename_);
        WriteSegments(stream);
        WriteLinks(stream);
        WritePaths(stream);
    }
};

struct ExtendedContigIdT {
    string full_id_;
    string short_id_;

    ExtendedContigIdT(): full_id_(""), short_id_("") {}

    ExtendedContigIdT(string full_id, string short_id): full_id_(full_id), short_id_(short_id) {}
};

template <class Graph>
void MakeContigIdMap(const Graph& graph, map<EdgeId, ExtendedContigIdT>& ids, const ConnectedComponentCounter &cc_counter_, string prefix) {
    int counter = 0;
    for (auto it = graph.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        if (ids.count(e) == 0) {
            string id;
            if (cfg::get().pd) {
                size_t c_id = cc_counter_.GetComponent(e);
                id = io::MakeContigComponentId(++counter, graph.length(e) + graph.k(), graph.coverage(e), c_id, prefix);
            }
            else
                id = io::MakeContigId(++counter, graph.length(e) + graph.k(), graph.coverage(e), prefix);
            ids[e] = ExtendedContigIdT(id, std::to_string(counter) + "+");
            if (e != graph.conjugate(e))
                ids[graph.conjugate(e)] =  ExtendedContigIdT(id + "'", std::to_string(counter) + "-");
        }
    }
}

//TODO refactor decently
template<class Graph>
class ContigPrinter {
private:
    const Graph &graph_;

    template<class sequence_stream>
    void ReportEdge(sequence_stream& oss, const string &s, double cov) {
        oss << s;
        oss << cov;
    }

    void ReportEdge(io::osequencestream_for_fastg& oss,
            const string& sequence,
            const string& id,
            const set<string>& nex_ids) {
        oss.set_header(id);
        oss << nex_ids;
        oss << sequence;
    }

public:
    ContigPrinter(const Graph &graph) : graph_(graph) {
    }

    template<class sequence_stream>
    void PrintContigs(sequence_stream &os) {
        for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            ReportEdge(os, graph_.EdgeNucls(e).str(), graph_.coverage(e));
        }
    }

    template<class sequence_stream>
    void PrintContigsFASTG(sequence_stream &os, const ConnectedComponentCounter & cc_counter) {
        map<EdgeId, ExtendedContigIdT> ids;
        MakeContigIdMap(graph_, ids, cc_counter, "EDGE");
        for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            set<string> next;
            for (EdgeId next_e : graph_.OutgoingEdges(graph_.EdgeEnd(e))) {
                next.insert(ids[next_e].full_id_);
            }
            ReportEdge(os, graph_.EdgeNucls(e).str(), ids[e].full_id_, next);
            if (e != graph_.conjugate(e)) {
                set<string> next_conj;
                for (EdgeId next_e : graph_.OutgoingEdges(graph_.EdgeEnd(graph_.conjugate(e)))) {
                    next_conj.insert(ids[next_e].full_id_);
                }
                ReportEdge(os, graph_.EdgeNucls(graph_.conjugate(e)).str(), ids[graph_.conjugate(e)].full_id_, next_conj);
            }
        }
    }
};

inline void OutputContigs(ConjugateDeBruijnGraph &g,
                          const string &contigs_output_filename) {
    INFO("Outputting contigs to " << contigs_output_filename << ".fasta");
    io::osequencestream_cov oss(contigs_output_filename + ".fasta");

    ContigPrinter<ConjugateDeBruijnGraph>(g).PrintContigs(oss);
}

inline void OutputContigsToGFA(ConjugateDeBruijnGraph &g, path_extend::PathContainer &paths, const string &contigs_output_filename) {
    INFO("Outputting graph to " << contigs_output_filename << ".gfa");
    GFAWriter<ConjugateDeBruijnGraph> writer(g, paths, contigs_output_filename + ".gfa");
    writer.Write();
}

inline void OutputContigsToFASTG(ConjugateDeBruijnGraph& g,
                   const string& contigs_output_filename, const ConnectedComponentCounter & cc_counter) {
    INFO("Outputting graph to " << contigs_output_filename << ".fastg");
    io::osequencestream_for_fastg ossfg(contigs_output_filename + ".fastg");
    ContigPrinter<ConjugateDeBruijnGraph>(g).PrintContigsFASTG(ossfg, cc_counter);
}

}
