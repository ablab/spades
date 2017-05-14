//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_support.hpp"

namespace path_extend {

template<class Graph>
using EdgeNamingF = std::function<std::string (const Graph&, EdgeId)>;

template<class Graph>
EdgeNamingF<Graph> IdNamingF(const string &prefix = "") {
    return [=](const Graph &g, EdgeId e) {
        return io::MakeContigId(g.int_id(e), prefix);
    };
}

template<class Graph>
EdgeNamingF<Graph> BasicNamingF(const string &prefix = "EDGE") {
    return [=](const Graph &g, EdgeId e) {
        return io::MakeContigId(g.int_id(e), g.length(e) + g.k(), g.coverage(e), prefix);
    };
}

template<class Graph>
EdgeNamingF<Graph> PlasmidNamingF(const ConnectedComponentCounter &cc_counter,
                           const string &prefix = "EDGE") {
    return [=, &cc_counter](const Graph &g, EdgeId e) {
        return io::MakeContigComponentId(g.int_id(e),
                                         g.length(e) + g.k(),
                                         g.coverage(e),
                                         cc_counter.GetComponent(e),
                                         prefix);
    };
}

template<class Graph>
class CanonicalEdgeHelper {
    const Graph &g_;
    const EdgeNamingF<Graph> naming_f_;
    const string pos_orient_;
    const string neg_orient_;
public:

    CanonicalEdgeHelper(const Graph &g,
                        EdgeNamingF<Graph> naming_f = IdNamingF<Graph>(),
                        const string& pos_orient = "+",
                        const string& neg_orient = "-") :
            g_(g), naming_f_(naming_f),
            pos_orient_(pos_orient), neg_orient_(neg_orient) {
    }

    bool IsCanonical(EdgeId e) const {
        return e <= g_.conjugate(e);
    }

    EdgeId Canonical(EdgeId e) const {
        return IsCanonical(e) ? e : g_.conjugate(e);
    }

    std::string GetOrientation(EdgeId e) const {
        return IsCanonical(e) ? pos_orient_ : neg_orient_;
    }

    std::string EdgeString(EdgeId e,
                           const std::string &delim = "") const {
        return naming_f_(g_, Canonical(e)) + delim + GetOrientation(e);
    }
};

template<class Graph>
class FastgWriter {
    typedef typename Graph::EdgeId EdgeId;

    const Graph &graph_;
    CanonicalEdgeHelper<Graph> short_namer_;
    CanonicalEdgeHelper<Graph> extended_namer_;

    string ToPathString(const BidirectionalPath &path) const {
        if (path.Empty())
            return "";
        string res = short_namer_.EdgeString(path.Front());
        for (size_t i = 1; i < path.Size(); ++i) {
            if (graph_.EdgeEnd(path[i - 1]) != graph_.EdgeStart(path[i]) || path.GapAt(i).gap > 0) {
                res += ";\n" + short_namer_.EdgeString(path[i]);
            } else {
                res += "," + short_namer_.EdgeString(path[i]);
            }
        }
        return res;
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

    FastgWriter(const Graph &graph,
                EdgeNamingF<Graph> edge_naming_f = BasicNamingF<Graph>())
            : graph_(graph),
              short_namer_(graph_),
              extended_namer_(graph_, edge_naming_f, "", "'") {
    }

    void WriteSegmentsAndLinks(const string &fn) {
        io::osequencestream_for_fastg os(fn);
        for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            set<string> next;
            for (EdgeId next_e : graph_.OutgoingEdges(graph_.EdgeEnd(e))) {
                next.insert(extended_namer_.EdgeString(next_e));
            }
            ReportEdge(os, graph_.EdgeNucls(e).str(), extended_namer_.EdgeString(e), next);
            if (e != graph_.conjugate(e)) {
                set<string> next_conj;
                for (EdgeId next_e : graph_.OutgoingEdges(graph_.EdgeEnd(graph_.conjugate(e)))) {
                    next_conj.insert(extended_namer_.EdgeString(next_e));
                }
                ReportEdge(os, graph_.EdgeNucls(graph_.conjugate(e)).str(),
                           extended_namer_.EdgeString(graph_.conjugate(e)), next_conj);
            }
        }
    }

    void WritePaths(const ScaffoldStorage &scaffold_storage, const string &fn) const {
        std::ofstream os(fn);
        for (const auto& scaffold_info : scaffold_storage) {
            os << scaffold_info.name << endl;
            os << ToPathString(*scaffold_info.path) << endl;
            os << scaffold_info.name << "'" << endl;
            os << ToPathString(*scaffold_info.path->GetConjPath()) << endl;
        }
    }

};

template<class Graph>
class GFAWriter {
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    CanonicalEdgeHelper<Graph> edge_namer_;
    std::ostream &os_;

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
        os_ << "L\t" << edge_namer_.EdgeString(e1, "\t") << "\t"
            << edge_namer_.EdgeString(e2, "\t") << "\t"
            << overlap_size << "M" << std::endl;
    }

    void WriteLinks() {
        for (VertexId v : graph_)
            for (auto inc_edge : graph_.IncomingEdges(v))
                for (auto out_edge : graph_.OutgoingEdges(v))
                    WriteLink(inc_edge, out_edge, graph_.k());
    }

    void WritePath(const std::string& name, size_t segment_id, const vector<std::string> &edge_strs) {
        os_ << "P" << "\t" ;
        os_ << name << "_" << segment_id << "\t";
        std::string delimeter = "";
        for (const auto& e : edge_strs) {
            os_ << delimeter << e;
            delimeter = ",";
        }
        os_ << "\t*\n";
//        delimeter = "";
//        for (size_t i = 0; i < edge_strs.size() - 1; ++i) {
//            os_ << delimeter << "*";
//            delimeter = ",";
//        }
//        os_ << "\n";
    }

public:
    GFAWriter(const Graph &graph, std::ostream &os,
              EdgeNamingF<Graph> naming_f = BasicNamingF<Graph>())
            : graph_(graph),
              os_(os),
              edge_namer_(graph_, naming_f) {
    }

    void WriteSegmentsAndLinks() {
        WriteSegments();
        WriteLinks();
    }

    void WritePaths(const ScaffoldStorage &scaffold_storage) {
        for (const auto& scaffold_info : scaffold_storage) {
            const path_extend::BidirectionalPath &p = *scaffold_info.path;
            if (p.Size() == 0) {
                continue;
            }
            std::vector<std::string> segmented_path;
            //size_t id = p.GetId();
            size_t segment_id = 1;
            for (size_t i = 0; i < p.Size() - 1; ++i) {
                EdgeId e = p[i];
                segmented_path.push_back(edge_namer_.EdgeString(e));
                if (graph_.EdgeEnd(e) != graph_.EdgeStart(p[i+1]) || p.GapAt(i+1).gap > 0) {
                    WritePath(scaffold_info.name, segment_id, segmented_path);
                    segment_id++;
                    segmented_path.clear();
                }
            }

            segmented_path.push_back(edge_namer_.EdgeString(p.Back()));
            WritePath(scaffold_info.name, segment_id, segmented_path);
        }
    }

};

void WriteScaffolds(const ScaffoldStorage &scaffold_storage, const string &fn);

typedef std::function<void (const ScaffoldStorage&)> PathsWriterT;

class ContigWriter {
    const Graph& g_;
    shared_ptr<ContigNameGenerator> name_generator_;

public:

    static PathsWriterT BasicFastaWriter(const string &fn) {
        return [=](const ScaffoldStorage& scaffold_storage) {
            WriteScaffolds(scaffold_storage, fn);
        };
    }

    ContigWriter(const Graph& g,
                 shared_ptr<ContigNameGenerator> name_generator) :
            g_(g),
            name_generator_(name_generator) {
    }

    void OutputPaths(const PathContainer &paths, const vector<PathsWriterT>& writers) const;

    void OutputPaths(const PathContainer &paths, PathsWriterT writer) const {
        OutputPaths(paths, vector<PathsWriterT>{writer});
    }

    void OutputPaths(const PathContainer &paths, const string &fn) const {
        OutputPaths(paths, BasicFastaWriter(fn));
    }

private:
    DECL_LOGGER("PathExtendIO")
};

}
