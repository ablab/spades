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

//This class corrects mismatches or masks repeat differences or other such things with the sequence of an edge
template<class Graph>
class ContigCorrector {
private:
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
protected:
    const Graph &graph() const {
        return graph_;
    }

public:
    ContigCorrector(const Graph &graph) : graph_(graph) {
    }

    virtual string correct(EdgeId e) = 0;

    virtual ~ContigCorrector() {
    }
};

template<class Graph>
class DefaultContigCorrector : public ContigCorrector<Graph> {
private:
    typedef typename Graph::EdgeId EdgeId;
public:
    DefaultContigCorrector(const Graph &graph) : ContigCorrector<Graph>(graph) {
    }

    string correct(EdgeId e) {
        return this->graph().EdgeNucls(e).str();
    }
};


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
private:
    std::ostream &ostream_;

public:

    GFAPathWriter(std::ostream &stream)
    : ostream_(stream)  {
    }

    void Write(const PathSegmentSequence &path_segment_sequence) {
        ostream_ << "P" << "\t" ;
        ostream_ << path_segment_sequence.path_id_ << "_" << path_segment_sequence.segment_number_ << "\t";
        std::string delimeter = "";
        for (size_t i = 0; i < path_segment_sequence.segment_sequence_.size() - 1; ++i) {
            ostream_ << delimeter << path_segment_sequence.segment_sequence_[i];
            delimeter = ",";
        }
        ostream_ << "\t";
        std::string delimeter2 = "";
        for (size_t i = 0; i < path_segment_sequence.segment_sequence_.size() - 1; ++i) {
                ostream_ << delimeter2 << "*";
                delimeter2 = ",";
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
    std::set<EdgeId> set_of_authentic_edges_;

    bool IsCanonical(EdgeId e) const {
        if (e <= graph_.conjugate(e)) {
            return true;
        } else {
            return false;
        }
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
                size_t segment_first = IsCanonical(inc_edge) ? inc_edge.int_id() : graph_.conjugate(inc_edge).int_id();
                for (auto out_edge : graph_.OutgoingEdges(*it)) {
                    size_t segment_second = IsCanonical(out_edge) ? out_edge.int_id() : graph_.conjugate(out_edge).int_id();
                    std::string orientation_second = GetOrientation(out_edge);
                    link_writer.Write(segment_first, orientation_first, segment_second, orientation_second);
                }
            }
        }
    }

    void UpdateSegmentedPath(PathSegmentSequence &segmented_path, EdgeId e) {
        std::string segment_id = IsCanonical(e) ? ToString(e.int_id()) : ToString(graph_.conjugate(e).int_id());
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

//This class uses corrected sequences to construct contig (just return as is, find unipath, trim contig)
template<class Graph>
class ContigConstructor {
private:
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    ContigCorrector<Graph> &corrector_;
protected:
    string correct(EdgeId e) {
        return corrector_.correct(e);
    }

    const Graph &graph() const {
        return graph_;
    }

public:

    ContigConstructor(const Graph &graph, ContigCorrector<Graph> &corrector) : graph_(graph), corrector_(corrector) {
    }

    virtual pair<string, double> construct(EdgeId e) = 0;

    virtual ~ContigConstructor(){
    }
};

template<class Graph>
class DefaultContigConstructor : public ContigConstructor<Graph> {
private:
    typedef typename Graph::EdgeId EdgeId;
public:

    DefaultContigConstructor(const Graph &graph, ContigCorrector<Graph> &corrector) : ContigConstructor<Graph>(graph, corrector) {
    }

    pair<string, double> construct(EdgeId e) {
        return make_pair(this->correct(e), this->graph().coverage(e));
    }
};

template<class Graph>
vector<typename Graph::EdgeId> Unipath(const Graph& g, typename Graph::EdgeId e) {
    omnigraph::UniquePathFinder<Graph> unipath_finder(g);
    vector<typename Graph::EdgeId> answer = unipath_finder.UniquePathBackward(e);
    const vector<typename Graph::EdgeId>& forward = unipath_finder.UniquePathForward(e);
    for (size_t i = 1; i < forward.size(); ++i) {
        answer.push_back(forward[i]);
    }
    return answer;
}

template<class Graph>
class UnipathConstructor : public ContigConstructor<Graph> {
private:
    typedef typename Graph::EdgeId EdgeId;



    string MergeOverlappingSequences(std::vector<string>& ss, size_t overlap) {
        if (ss.empty()) {
            return "";
        }
        stringstream result;
        result << ss.front().substr(0, overlap);
//        prev_end = ss.front().substr(0, overlap);
        for (auto it = ss.begin(); it != ss.end(); ++it) {
//            VERIFY(prev_end == it->substr(0, overlap));
            result << it->substr(overlap);
//            prev_end = it->substr(it->size() - overlap);
        }
        return result.str();
    }


    string MergeSequences(const Graph& g,
            const vector<typename Graph::EdgeId>& continuous_path) {
        vector<string> path_sequences;
        for (size_t i = 0; i < continuous_path.size(); ++i) {
            if(i > 0)
                VERIFY(
                    g.EdgeEnd(continuous_path[i - 1])
                            == g.EdgeStart(continuous_path[i]));
            path_sequences.push_back(this->correct(continuous_path[i]));
        }
        return MergeOverlappingSequences(path_sequences, g.k());
    }

public:

    UnipathConstructor(const Graph &graph, ContigCorrector<Graph> &corrector) : ContigConstructor<Graph>(graph, corrector) {
    }

    pair<string, double> construct(EdgeId e) {
        vector<EdgeId> unipath = Unipath(this->graph(), e);
        return make_pair(MergeSequences(this->graph(), unipath), stats::AvgCoverage(this->graph(), unipath));
    }
};

template<class Graph>
class CuttingContigConstructor : public ContigConstructor<Graph> {
private:
    typedef typename Graph::EdgeId EdgeId;

    bool ShouldCut(VertexId v) const {
        const Graph &g = this->graph();
        vector<EdgeId> edges;
        push_back_all(edges, g.OutgoingEdges(v));
        if(edges.size() == 0)
            return false;
        for(size_t i = 1; i < edges.size(); i++) {
            if(g.EdgeNucls(edges[i])[g.k()] != g.EdgeNucls(edges[0])[g.k()])
                return false;
        }
        edges.clear();
        push_back_all(edges, g.IncomingEdges(v));
        for(size_t i = 0; i < edges.size(); i++)
            for(size_t j = i + 1; j < edges.size(); j++) {
                if(g.EdgeNucls(edges[i])[g.length(edges[i]) - 1] != g.EdgeNucls(edges[j])[g.length(edges[j]) - 1])
                    return true;
            }
        return false;
    }

public:

    CuttingContigConstructor(const Graph &graph, ContigCorrector<Graph> &corrector) : ContigConstructor<Graph>(graph, corrector) {
    }

    pair<string, double> construct(EdgeId e) {
        string result = this->correct(e);
        if(result.size() > this->graph().k() && ShouldCut(this->graph().EdgeEnd(e))) {
            result = result.substr(0, result.size() - this->graph().k());
        }
        if(result.size() > this->graph().k() && ShouldCut(this->graph().conjugate(this->graph().EdgeStart(e)))) {
            result = result.substr(this->graph().k(), result.size());
        }
        return make_pair(result, this->graph().coverage(e));
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
            ids[e] = ExtendedContigIdT(id, ToString(counter) + "+");
            if (e != graph.conjugate(e))
                ids[graph.conjugate(e)] =  ExtendedContigIdT(id + "'", ToString(counter) + "-");
        }
    }
}

template<class Graph>
class ContigPrinter {
private:
    const Graph &graph_;
    ContigConstructor<Graph> &constructor_;
    template<class sequence_stream>
    void ReportEdge(sequence_stream& oss
            , const pair<string, double> sequence_data) {
        oss << sequence_data.second;
        oss << sequence_data.first;
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
    ContigPrinter(const Graph &graph, ContigConstructor<Graph> &constructor) : graph_(graph), constructor_(constructor) {
    }

    template<class sequence_stream>
    void PrintContigs(sequence_stream &os) {
        for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            ReportEdge<sequence_stream>(os, constructor_.construct(*it));
        }
    }

    template<class sequence_stream>
    void PrintContigsFASTG(sequence_stream &os, const ConnectedComponentCounter & cc_counter) {
        map<EdgeId, ExtendedContigIdT> ids;
        MakeContigIdMap(graph_, ids, cc_counter, "EDGE");
        for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            set<string> next;
            VertexId v = graph_.EdgeEnd(*it);
            auto edges = graph_.OutgoingEdges(v);
            for (auto next_it = edges.begin(); next_it != edges.end(); ++next_it) {
                next.insert(ids[*next_it].full_id_);
            }
            ReportEdge(os, constructor_.construct(*it).first, ids[*it].full_id_, next);
            if (*it != graph_.conjugate(*it))
            {
                set<string> next_conj;
                v = graph_.EdgeEnd(graph_.conjugate(*it));
                edges = graph_.OutgoingEdges(v);
                for (auto next_it = edges.begin(); next_it != edges.end(); ++next_it) {
                    next_conj.insert(ids[*next_it].full_id_);
                }
                ReportEdge(os, constructor_.construct(graph_.conjugate(*it)).first, ids[graph_.conjugate(*it)].full_id_, next_conj);               
            }
        }
    }
};

template<class Graph>
bool PossibleECSimpleCheck(const Graph& g
        , typename Graph::EdgeId e) {
    return g.OutgoingEdgeCount(g.EdgeStart(e)) > 1 && g.IncomingEdgeCount(g.EdgeEnd(e)) > 1;
}

template<class Graph>
void ReportEdge(io::osequencestream_cov& oss
        , const Graph& g
        , typename Graph::EdgeId e
        , bool output_unipath = false
        , size_t solid_edge_length_bound = 0) {
    typedef typename Graph::EdgeId EdgeId;
    if (!output_unipath || (PossibleECSimpleCheck(g, e) && g.length(e) <= solid_edge_length_bound)) {
        TRACE("Outputting edge " << g.str(e) << " as single edge");
        oss << g.coverage(e);
        oss << g.EdgeNucls(e);
    } else {
        TRACE("Outputting edge " << g.str(e) << " as part of unipath");
        vector<EdgeId> unipath = Unipath(g, e);
        TRACE("Unipath is " << g.str(unipath));
        oss << stats::AvgCoverage(g, unipath);
        TRACE("Merged sequence is of length " << MergeSequences(g, unipath).size());
        oss << MergeSequences(g, unipath);
    }
}

inline void OutputContigs(ConjugateDeBruijnGraph &g, const string &contigs_output_filename, bool output_unipath) {
    INFO("Outputting contigs to " << contigs_output_filename << ".fasta");
    DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(g);
    io::osequencestream_cov oss(contigs_output_filename + ".fasta");

    if(!output_unipath) {
        DefaultContigConstructor<ConjugateDeBruijnGraph> constructor(g, corrector);

        ContigPrinter<ConjugateDeBruijnGraph>(g, constructor).PrintContigs(oss);
    } else {
        UnipathConstructor<ConjugateDeBruijnGraph> constructor(g, corrector);
        ContigPrinter<ConjugateDeBruijnGraph>(g, constructor).PrintContigs(oss);
    }

//    {
//        osequencestream_cov oss(contigs_output_filename);
//        set<ConjugateDeBruijnGraph::EdgeId> edges;
//        for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
//            if (edges.count(*it) == 0) {
//                ReportEdge(oss, g, *it, output_unipath, solid_edge_length_bound + ".oppa.fasta");
//                edges.insert(g.conjugate(*it));
//            }
//            //        oss << g.EdgeNucls(*it);
//        }
//        DEBUG("Contigs written");
//    }
//    if(!output_unipath) {
//        OutputContigs(g, contigs_output_filename + ".2.fasta", true, solid_edge_length_bound);
//    }
}

inline void OutputContigsToGFA(ConjugateDeBruijnGraph &g, path_extend::PathContainer &paths, const string &contigs_output_filename) {
    INFO("Outputting graph to " << contigs_output_filename << ".gfa");
    GFAWriter<ConjugateDeBruijnGraph> writer(g, paths, contigs_output_filename + ".gfa");
    writer.Write();
}


inline void OutputContigsToFASTG(ConjugateDeBruijnGraph& g,
                   const string& contigs_output_filename, const ConnectedComponentCounter & cc_counter) {

    INFO("Outputting graph to " << contigs_output_filename << ".fastg");
    DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(g);
    DefaultContigConstructor<ConjugateDeBruijnGraph> constructor(g, corrector);
    io::osequencestream_for_fastg ossfg(contigs_output_filename + ".fastg");
    ContigPrinter<ConjugateDeBruijnGraph>(g, constructor).PrintContigsFASTG(ossfg, cc_counter);
}




inline bool ShouldCut(ConjugateDeBruijnGraph& g, VertexId v) {
    vector<EdgeId> edges;
    push_back_all(edges, g.OutgoingEdges(v));

    if(edges.size() == 0)
        return false;
    for(size_t i = 1; i < edges.size(); i++) {
        if(g.EdgeNucls(edges[i])[g.k()] != g.EdgeNucls(edges[0])[g.k()])
            return false;
    }
    edges.clear();
    push_back_all(edges, g.IncomingEdges(v));
    for(size_t i = 0; i < edges.size(); i++)
        for(size_t j = i + 1; j < edges.size(); j++) {
            if(g.EdgeNucls(edges[i])[g.length(edges[i]) - 1] != g.EdgeNucls(edges[j])[g.length(edges[j]) - 1])
                return true;
        }
    return false;
}

inline void OutputCutContigs(ConjugateDeBruijnGraph& g,
        const string& contigs_output_filename,
        bool /*output_unipath*/ = false,
        size_t /*solid_edge_length_bound*/ = 0) {
    INFO("Outputting contigs to " << contigs_output_filename);
    DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(g);
    io::osequencestream_cov oss(contigs_output_filename);
    CuttingContigConstructor<ConjugateDeBruijnGraph> constructor(g, corrector);

//    osequencestream_cov oss(contigs_output_filename);
//    set<ConjugateDeBruijnGraph::EdgeId> edges;
//    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
//        EdgeId e = *it;
//        cout << g.length(e) << endl;
//        if (edges.count(e) == 0) {
//            Sequence s = g.EdgeNucls(e);
//            cout << s.size() << endl;
//            cout << "oppa " << ShouldCut(g, g.EdgeEnd(e)) << endl;
//            if(s.size() > g.k() && ShouldCut(g, g.EdgeEnd(e))) {
//                s = s.Subseq(0, s.size() - g.k());
//                cout << s.size() << endl;
//            }
//            cout << "oppa1 " << ShouldCut(g, g.conjugate(g.EdgeStart(e))) << endl;
//            if(s.size() > g.k() && ShouldCut(g, g.conjugate(g.EdgeStart(e)))) {
//                s = s.Subseq(g.k(), s.size());
//                cout << s.size() << endl;
//            }
//            oss << g.coverage(e);
//            oss << s;
//            edges.insert(g.conjugate(*it));
//        }
//        //        oss << g.EdgeNucls(*it);
//    }
}

inline void OutputSingleFileContigs(ConjugateDeBruijnGraph& g,
        const string& contigs_output_dir) {
    INFO("Outputting contigs to " << contigs_output_dir);
    int n = 0;
    make_dir(contigs_output_dir);
    char n_str[20];
    set<ConjugateDeBruijnGraph::EdgeId> edges;
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        if (edges.count(*it) == 0) {
            sprintf(n_str, "%d.fa", n);
            edges.insert(g.conjugate(*it));
            io::osequencestream oss(contigs_output_dir + n_str);
            oss << g.EdgeNucls(*it);
            n++;
        }
    }
    DEBUG("SingleFileContigs(Conjugate) written");
}

}
