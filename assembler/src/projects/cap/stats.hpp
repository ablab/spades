//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/graph_support/component_filters.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "assembly_graph/components/splitters.hpp"
#include "utils.hpp"
#include "utils/simple_tools.hpp"
#include "comparison_utils.hpp"
#include "assembly_graph/graph_support/basic_graph_stats.hpp"
#include "coloring.hpp"
#include "visualization/visualization_utils.hpp"

namespace cap {

template<class Graph>
class Component {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    vector<size_t> edge_lengths_;
    vector<VertexId> component_;
public:
    Component(const Graph &g, const vector<VertexId> &component) :
            component_(component) {
        for (size_t i = 0; i < component.size(); i++)
            for (size_t j = 0; j < component.size(); j++) {
                vector<EdgeId> edges = g.GetEdgesBetween(component[i],
                        component[j]);
                for (auto it = edges.begin(); it != edges.end(); ++it) {
                    edge_lengths_.push_back(g.length(*it));
                }
            }
        std::sort(edge_lengths_.rbegin(), edge_lengths_.rend());
    }

    bool operator<(const Component<Graph> &that) const {
        size_t i = 0;
        while (i < this->edge_lengths_.size() && i < that.edge_lengths_.size()
                && this->edge_lengths_[i] == that.edge_lengths_[i])
            i++;
        if (i == that.edge_lengths_.size())
            return false;
        if (i == this->edge_lengths_.size())
            return true;
        return this->edge_lengths_[i] < that.edge_lengths_[i];
    }

    bool operator==(const Component<Graph> &that) const {
        if (this->edge_lengths_.size() != that.edge_lengths_.size())
            return false;
        for (size_t i = 0; i < this->edge_lengths_.size(); i++)
            if (this->edge_lengths_[i] != that.edge_lengths_[i])
                return false;
        return true;
    }

    const vector<size_t> &edge_lengths() const {
        return edge_lengths_;
    }

};

template<class Stream, class Graph>
Stream &operator<<(Stream &stream, const Component<Graph> &component) {
    const vector<size_t> &lengths = component.edge_lengths();
    for (size_t i = 0; i < lengths.size(); i++) {
        stream << lengths[i] << " ";
    }
    return stream;
}

template<class Graph>
class ComponentClassifier {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

public:
    enum component_type {
        any = 0,
        error,
        single_red,
        single_blue,
        simple_bulge,
        tip,
        simple_misassembly,
        monochrome,
        complex_misassembly,
        size
    };

    static string info_printer_pos_name(size_t t) {
        const string names[] = { "all", "error", "single_red", "single_blue",
                "simple_bulge", "tip", "simple_misassembly", "monochrome",
                "complex_misassembly", "size" };
        return names[t];
    }

private:
    const Graph &g_;
    const ColorHandler<Graph> &coloring_;
    const size_t bulge_length_;

public:
    ComponentClassifier(const Graph &g, const ColorHandler<Graph> &coloring,
            size_t bulge_length) :
            g_(g), coloring_(coloring), bulge_length_(bulge_length) {
    }

    ComponentClassifier(const Graph &g, const ColorHandler<Graph> &coloring) :
            g_(g), coloring_(coloring), bulge_length_(g_.k() * 1000000) {
    }

    TColorSet GetColour(EdgeId edge) const {
        return coloring_.Color(edge);
    }

    bool CheckSimpleMisassembly(const vector<VertexId> &component) const {
        if (component.size() == 4) {
            for (size_t i = 0; i < 4; i++)
                for (size_t j = i + 1; j < 4; j++) {
                    vector<VertexId> sources;
                    sources.push_back(component[i]);
                    sources.push_back(component[j]);
                    vector<VertexId> sinks;
                    for (size_t k = 0; k < 4; k++) {
                        if (k != i && k != j)
                            sinks.push_back(component[k]);
                    }
                    if (CheckSimpleMisassembly(sources, sinks)) {
                        return true;
                    }
                }
        }
        return false;
    }

    bool CheckSimpleMisassembly(const vector<VertexId>& sources,
            const vector<VertexId>& sinks) const {
        if (sources.size() != 2 || sinks.size() != 2)
            return false;
        for (size_t i = 0; i < sources.size(); i++)
            for (size_t j = 0; j < sinks.size(); j++) {
                if (g_.GetEdgesBetween(sources[i], sinks[j]).size() != 1) {
                    return false;
                }
            }
        for (size_t i = 0; i < 2; i++) {
            if (g_.GetEdgesBetween(sources[i], sources[1 - i]).size() != 0)
                return false;
            if (g_.GetEdgesBetween(sinks[i], sinks[1 - i]).size() != 0)
                return false;
        }
        for (size_t i = 0; i < 2; i++) {
            if (GetColour(g_.GetEdgesBetween(sources[i], sinks[0])[0])
                    == GetColour(g_.GetEdgesBetween(sources[i], sinks[1])[0]))
                return false;
            if (GetColour(g_.GetEdgesBetween(sources[0], sinks[i])[0])
                    == GetColour(g_.GetEdgesBetween(sources[1], sinks[i])[0]))
                return false;
        }
        return true;
    }

    bool CheckIsolated(TColorSet colour,
            const vector<VertexId> &component) const {
        if (component.size() != 2)
            return false;
        vector<EdgeId> edges01 = g_.GetEdgesBetween(component[0], component[1]);
        vector<EdgeId> edges10 = g_.GetEdgesBetween(component[1], component[0]);
        vector<EdgeId> edges;
        edges.insert(edges.end(), edges01.begin(), edges01.end());
        edges.insert(edges.end(), edges10.begin(), edges10.end());
        if (edges.size() != 1) {
            return false;
        }
        return GetColour(edges[0]) == colour;
    }

    bool CheckBulge(const vector<VertexId> &component) const {
        if (component.size() != 2)
            return false;
        vector<EdgeId> edges01 = g_.GetEdgesBetween(component[0], component[1]);
        vector<EdgeId> edges10 = g_.GetEdgesBetween(component[1], component[0]);
        vector<EdgeId> edges;
        edges.insert(edges.end(), edges01.begin(), edges01.end());
        edges.insert(edges.end(), edges10.begin(), edges10.end());
        return (edges01.size() == 0 || edges10.size() == 0) && edges.size() == 2
                && g_.length(edges[0]) < bulge_length_
                && g_.length(edges[1]) < bulge_length_;
    }

    size_t EdgeNumber(const vector<VertexId> &component) const {
        size_t result = 0;
        for (size_t i = 0; i < component.size(); i++)
            for (size_t j = 0; j < component.size(); j++) {
                result += g_.GetEdgesBetween(component[i], component[j]).size();
            }
        return result;
    }

    bool Connected(VertexId v1, VertexId v2) const {
        return g_.GetEdgesBetween(v1, v2).size() > 0;
    }

    bool CheckTip(const vector<VertexId> &component) const {
        if (component.size() != 3)
            return false;
        if (EdgeNumber(component) != 2)
            return false;
        for (size_t i = 0; i < 3; i++) {
            if (CheckFork(component[i], component[(i + 1) % 3],
                    component[(i + 2) % 3]))
                return true;
        }
        return false;
    }

    bool CheckFork(VertexId base, VertexId tip1, VertexId tip2) const {
        return (Connected(base, tip1) && Connected(base, tip2))
                || (Connected(tip1, base) && Connected(tip2, base));
    }

    bool CheckMonochrome(const vector<VertexId> &component) const {
        set<TColorSet> colours;
        for (size_t i = 0; i < component.size(); i++)
            for (size_t j = 0; j < component.size(); j++) {
                vector<EdgeId> edges = g_.GetEdgesBetween(component[i],
                        component[j]);
                for (auto it = edges.begin(); it != edges.end(); ++it) {
                    colours.insert(GetColour(*it));
                }
            }
        return colours.size() == 1;
    }

    component_type GetComponentType(const vector<VertexId> &component) const {
        if (component.size() < 2)
            return component_type::error;
        if (component.size() == 2) {
            if (CheckIsolated(kRedColorSet, component))
                return component_type::single_red;
            if (CheckIsolated(kBlueColorSet, component))
                return component_type::single_blue;
            if (CheckBulge(component)) {
                return component_type::simple_bulge;
            }
            return component_type::complex_misassembly;
        }
        if (CheckMonochrome(component))
            return component_type::monochrome;
        if (CheckTip(component)) {
            return component_type::tip;
        }
        if (CheckSimpleMisassembly(component))
            return component_type::simple_misassembly;
        return component_type::complex_misassembly;
    }
};

template<class Graph>
class ComponentTypeFilter: public GraphComponentFilter<Graph> {
private:
    typedef GraphComponentFilter<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef typename ComponentClassifier<Graph>::component_type component_type;

    component_type to_draw_;
    ComponentClassifier<Graph> cc_;

public:
    ComponentTypeFilter(const Graph &g, component_type to_draw,
                        const ColorHandler<Graph> &coloring)
            : base(g),
              to_draw_(to_draw),
              cc_(g, coloring) {
    }

    /*virtual*/
    bool Check(const vector<VertexId>& component) const {
        return to_draw_ == component_type::any
                || cc_.GetComponentType(component) == to_draw_;
    }

private:
    DECL_LOGGER("ComponentTypeFilter")
    ;
};

template<class Graph>
class BreakPointsFilter: public GraphComponentFilter<Graph> {
    typedef GraphComponentFilter<Graph> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const ColorHandler<Graph> coloring_;
    size_t color_threshold_;
public:
    BreakPointsFilter(const Graph& graph, const ColorHandler<Graph>& coloring,
            size_t color_threshold) :
            base(graph), coloring_(coloring), color_threshold_(color_threshold) {

    }

    bool MultiColored(const GraphComponent<Graph>& component) const {
        set<TColorSet> colors;
        for (auto it = component.e_begin(); it != component.e_end(); ++it) {
            colors.insert(coloring_.Color(*it));
        }
        return colors.size() >= color_threshold_;
    }

    /*virtual*/
    //todo change to set or GraphComponent and add useful protected methods
    bool Check(const vector<VertexId> &component_veritces) const {
        GraphComponent<Graph> component(this->graph(),
                component_veritces.begin(), component_veritces.end());
        return component.v_size() > 2 && MultiColored(component);
    }

};

template<class Graph>
class BreakPointGraphStatistics: public GraphComponentFilter<Graph> {
private:
    typedef GraphComponentFilter<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef typename ComponentClassifier<Graph>::component_type component_type;

    const ColorHandler<Graph> &coloring_;
    ComponentClassifier<Graph> cc_;
    mutable vector<vector<Component<Graph>>> components_;
    mutable vector<size_t> total_red_;
    mutable vector<size_t> total_blue_;
    bool ready_;

    void UpdateTotalEdgeLength(component_type t, const vector<VertexId>& component) const {
        for(size_t i = 0; i < component.size(); i++) {
            for(size_t j = 0; j < component.size(); j++) {
                vector<EdgeId> edges = this->graph().GetEdgesBetween(component[i], component[j]);
                for(auto it = edges.begin(); it != edges.end(); ++it) {
                    if(coloring_.Color(*it) == kRedColorSet)
                    total_red_[t] += this->graph().length(*it);
                    if(coloring_.Color(*it) == kBlueColorSet)
                    total_blue_[t] += this->graph().length(*it);
                }
            }
        }
    }

    void UpdateStats(const vector<VertexId>& component) const {
        component_type t = cc_.GetComponentType(component);
        Component<Graph> c(this->graph(), component);
        components_[t].push_back(c);
        UpdateTotalEdgeLength(t, component);
    }

public:
    BreakPointGraphStatistics(const Graph &g, const ColorHandler<Graph> &coloring) :
    base(g), coloring_(coloring), cc_(g, coloring), components_(component_type::size), total_red_(component_type::size), total_blue_(component_type::size), ready_(
            false) {
    }

    /*virtual*/bool Check(const vector<VertexId>& component) const {
        UpdateStats(component);
        return false;
    }

    void CountStats() {
        visualization::graph_labeler::EmptyGraphLabeler<Graph> labeler;
        make_dir("assembly_compare");
        shared_ptr<GraphSplitter<Graph>> splitter = LongEdgesExclusiveSplitter<Graph>(this->graph(), 1000000000);
        WriteComponents(this->graph(), *splitter, *this,
                "assembly_compare/breakpoint_graph.dot",
                *visualization::ConstructColorer(coloring_), labeler);
        ready_ = true;
        for (size_t i = 0; i < component_type::size; ++i) {
            INFO("Number of components of type " << ComponentClassifier<Graph>::info_printer_pos_name(i) << " is " << GetComponentNumber(i));
            INFO("Total length of red edges in these components is " << total_red_[i]);
            INFO("Total length of blue edges in these components is " << total_blue_[i]);
        }
    }

    size_t GetComponentNumber(size_t t) const {
        return components_[t].size();
    }

    const vector<Component<Graph>> &GetComponents(size_t t) const {
        std::sort(components_[t].rbegin(), components_[t].rend());
        return components_[t];
    }

private:
    DECL_LOGGER("BreakPointGraphStatistics");
};

template<class Graph>
class BPGraphStatCounter {
private:
    typedef typename ComponentClassifier<Graph>::component_type component_type;
    typedef typename Graph::VertexId VertexId;
    const Graph &graph_;
    const ColorHandler<Graph> &coloring_;
    const string output_folder_;
public:
    BPGraphStatCounter(const Graph &g, const ColorHandler<Graph> &coloring,
            const string& output_folder) :
            graph_(g), coloring_(coloring), output_folder_(output_folder) {
    }

    void PrintComponents(component_type c_type,
            const visualization::graph_labeler::GraphLabeler<Graph>& labeler,
            bool create_subdir = true) const {
        string filename;
        if (create_subdir) {
            make_dir(output_folder_);
            string type_dir = output_folder_
                    + ComponentClassifier<Graph>::info_printer_pos_name(c_type)
                    + "/";
            make_dir(type_dir);
            string picture_dir = type_dir + "pictures/";
            make_dir(picture_dir);
            filename = picture_dir + "breakpoint_graph.dot";
        } else {
            filename = output_folder_
                    + ComponentClassifier<Graph>::info_printer_pos_name(c_type)
                    + ".dot";
        }
        shared_ptr<GraphSplitter<Graph>> splitter = LongEdgesExclusiveSplitter<Graph>(graph_, 1000000000);
        ComponentTypeFilter<Graph> stats(graph_, c_type, coloring_);

        WriteComponents(this->graph_, *splitter, stats, filename,
                *ConstructColorer(coloring_), labeler);
    }

    void PrintStats(const BreakPointGraphStatistics<Graph> &stats) const {
        make_dir(output_folder_);
        for (size_t t = 0; t < component_type::size; t++) {
            string type_dir = output_folder_
                    + ComponentClassifier<Graph>::info_printer_pos_name(t)
                    + "/";
            make_dir(type_dir);
            ofstream stream;
            stream.open((type_dir + "components.txt").c_str());
            const vector<Component<Graph>> &components = stats.GetComponents(t);
            for (auto it = components.begin(); it != components.end(); ++it) {
                stream << *it << endl;
            }
            stream.close();
        }
    }

    void CountStats(const visualization::graph_labeler::GraphLabeler<Graph>& labeler, bool detailed_output =
            true) const {
        make_dir(output_folder_);
        BreakPointGraphStatistics<Graph> stats(graph_, coloring_);
        stats.CountStats();
        if (detailed_output) {
            PrintStats(stats);
            PrintComponents(component_type::complex_misassembly, labeler);
            PrintComponents(component_type::monochrome, labeler);
            PrintComponents(component_type::tip, labeler);
            PrintComponents(component_type::simple_misassembly, labeler);
        }
        PrintComponents(component_type::any, labeler, detailed_output);
    }
};

template<class Graph>
class TrivialBreakpointFinder: public AbstractFilter<
        vector<typename Graph::VertexId>> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    struct bp_comp {
        bp_comp(const Graph& g, const ColorHandler<Graph>& coloring) :
                g_(g), coloring_(coloring) {
        }

        bool operator()(VertexId v1, VertexId v2) {
            return MaxRedBlueIncLength(v1) > MaxRedBlueIncLength(v2);
        }

        size_t MaxRedBlueIncLength(VertexId v) {
            vector<EdgeId> edges;
            push_back_all(edges, g_.IncomingEdges(v));
            push_back_all(edges, g_.OutgoingEdges(v));
            return MaxRedBlueLength(edges);
        }

    private:
        const Graph& g_;
        const ColorHandler<Graph>& coloring_;

        size_t MaxRedBlueLength(const vector<EdgeId> edges) {
            size_t max_length = 0;
            for (auto it = edges.begin(); it != edges.end(); ++it) {
                if (coloring_.Color(*it) == kBlueColorSet
                        || coloring_.Color(*it) == kRedColorSet) {
                    if (g_.length(*it) > max_length) {
                        max_length = g_.length(*it);
                    }
                }
            }
            VERIFY(max_length > 0);
            return max_length;
        }
    };

    const Graph& g_;
    const ColorHandler<Graph>& coloring_;
    const EdgesPositionHandler<Graph>& pos_;

    void ReportBreakpoint(VertexId v, const string& folder,
            const string& prefix) {
        TRACE("Vertex " << g_.str(v) << " identified as breakpoint");
        visualization::graph_labeler::LengthIdGraphLabeler<Graph> basic_labeler(g_);
        visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(g_, pos_);

        visualization::graph_labeler::CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
        VERIFY(g_.OutgoingEdgeCount(v) > 0);
        EdgeId e = g_.OutgoingEdges(v).front();
        GraphComponent<Graph> component = omnigraph::EdgeNeighborhood(g_, e);
        visualization::visualization_utils::WriteComponent(
                component,
                folder + prefix + ToString(g_.int_id(v)) + "_loc.dot",
                coloring_.ConstructColorer(component), labeler);
    }

    bool CheckEdges(const vector<EdgeId>& edges) const {
        std::set<TColorSet> colors;
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            colors.insert(coloring_.Color(*it));
        }
        return edges.size() == 2 && colors.count(kBlueColorSet) == 1
                && NotTips(edges);
    }

    bool IsTip(VertexId v) const {
        return g_.IncomingEdgeCount(v) + g_.OutgoingEdgeCount(v) == 1;
    }

    bool IsTip(EdgeId e) const {
        return (IsTip(g_.EdgeStart(e)) || IsTip(g_.EdgeEnd(e)))
                && g_.length(e) < 200;
    }

    bool NotTips(const vector<EdgeId>& edges) const {
        for (auto it = edges.begin(); it != edges.end(); ++it)
            if (IsTip(*it))
                return false;
        return true;
    }

public:
    TrivialBreakpointFinder(const Graph& g, const ColorHandler<Graph>& coloring,
            const EdgesPositionHandler<Graph>& pos) :
            g_(g), coloring_(coloring), pos_(pos) {
    }

    void FindBreakPoints(const string& folder) {
        vector<VertexId> breakpoints;
        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            if (coloring_.Color(*it) == kRedColorSet) {
                if (CheckEdges(g_.OutgoingEdges(g_.EdgeStart(*it))))
                    breakpoints.push_back(g_.EdgeStart(*it));
//                    ReportBreakpoint(g_.EdgeStart(*it));
                if (CheckEdges(g_.IncomingEdges(g_.EdgeEnd(*it))))
                    breakpoints.push_back(g_.EdgeEnd(*it));
//                    ReportBreakpoint(g_.EdgeEnd(*it));
            }
        }
        bp_comp comp(g_, coloring_);
        std::sort(breakpoints.begin(), breakpoints.end(), comp);
        for (size_t i = 0; i < breakpoints.size(); ++i) {
            ReportBreakpoint(
                    breakpoints[i],
                    folder,
                    ToString(i) + "_"
                            + ToString(comp.MaxRedBlueIncLength(breakpoints[i]))
                            + "_");
        }
    }

    virtual bool Check(const vector<typename Graph::VertexId> &vertices) const {
        GraphComponent<Graph> component(g_, vertices.begin(), vertices.end());
        for (auto it = component.e_begin(); it != component.e_end(); ++it) {
            if (coloring_.Color(*it) == kRedColorSet) {
                if (CheckEdges(g_.OutgoingEdges(g_.EdgeStart(*it)))
                        || CheckEdges(g_.IncomingEdges(g_.EdgeEnd(*it))))
                    return true;
            }
        }
        return false;
    }

private:
    DECL_LOGGER("TrivialBreakpointFinder")
    ;
};

template<class Graph>
class SimpleInDelAnalyzer {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph& g_;
    const ColorHandler<Graph>& coloring_;
    const EdgesPositionHandler<Graph>& edge_pos_;
    const vector<EdgeId> genome_path_;
    const TColorSet shortcut_color_;
    const string folder_;

    vector<EdgeId> TryFindPath(size_t pos, VertexId end,
            size_t edge_count_bound) {
        vector<EdgeId> answer;
        for (size_t i = 0;
                i + pos < genome_path_.size() && i < edge_count_bound; ++i) {
//            if ((coloring_.Color(genome_path_[pos + i]) & shortcut_color_) > 0) {
//                DEBUG("Came into edge of wrong color");
//                return vector<EdgeId>();
//            }
            answer.push_back(genome_path_[pos + i]);
            if (g_.EdgeEnd(genome_path_[pos + i]) == end) {
                return answer;
            }
        }
        DEBUG("Edge bound reached");
        return vector<EdgeId>();
    }

    map<TColorSet, size_t> ColorLengths(const vector<EdgeId>& edges) {
        map<TColorSet, size_t> answer;
        for (size_t i = 0; i < edges.size(); ++i) {
            answer[coloring_.Color(edges[i])] += g_.length(edges[i]);
        }
        return answer;
    }

    size_t VioletLengthOfGenomeUnique(const vector<EdgeId>& edges) {
        size_t answer = 0;
        for (size_t i = 0; i < edges.size(); ++i) {
            if (coloring_.Color(edges[i]) == kVioletColorSet
                    && std::count(genome_path_.begin(), genome_path_.end(), edges[i]) == 1) {
                answer += g_.length(edges[i]);
            }
        }
        return answer;
    }

    //genome pos exclusive
    size_t CumulativeGenomeLengthToPos(size_t pos) {
        size_t answer = 0;
        for (size_t i = 0; i < pos; ++i) {
            answer += g_.length(genome_path_[i]);
        }
        return answer;
    }

    pair<vector<EdgeId>, pair<size_t, size_t>> FindGenomePath(VertexId start,
            VertexId end, size_t edge_count_bound) {
        for (size_t i = 0; i < genome_path_.size(); ++i) {
            if (g_.EdgeStart(genome_path_[i]) == start) {
                vector<EdgeId> path = TryFindPath(i, end, edge_count_bound);
                if (!path.empty())
                    return make_pair(
                            path,
                            make_pair(
                                    CumulativeGenomeLengthToPos(i),
                                    CumulativeGenomeLengthToPos(
                                            i + path.size())));
            }
        }
        return make_pair(vector<EdgeId>(), make_pair(0, 0));
    }

    pair<string, pair<size_t, size_t>> ContigIdAndPositions(EdgeId e) {
        vector<EdgePosition> poss = edge_pos_.GetEdgePositions(e);
        VERIFY(!poss.empty());
        if (poss.size() > 1) {
            WARN("Something strange with assembly positions");
            return make_pair("", make_pair(0, 0));
        }
        EdgePosition pos = poss.front();
        return make_pair(pos.contigId, make_pair(pos.mr.initial_range.start_pos, pos.mr.initial_range.end_pos));
    }

    void WriteAltPath(EdgeId e, const vector<EdgeId>& genome_path) {
        visualization::graph_labeler::LengthIdGraphLabeler<Graph> basic_labeler(g_);
        visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(g_, edge_pos_);

        visualization::graph_labeler::CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);

        string alt_path_folder = folder_ + ToString(g_.int_id(e)) + "/";
        make_dir(alt_path_folder);
        WriteComponentsAlongPath(g_, labeler, alt_path_folder + "path.dot", /*split_length*/
                1000, /*vertex_number*/15, TrivialMappingPath(g_, genome_path),
                *ConstructBorderColorer(g_, coloring_));
    }

    void Process(EdgeId e, const vector<EdgeId>& genome_path,
            size_t genome_start, size_t genome_end) {
        DEBUG("Processing edge and genome path");
        const size_t mem_lim = 2 << 26;
        Sequence edge_nucls = g_.EdgeNucls(e);
        Sequence path_nucls = debruijn_graph::MergeSequences(g_, genome_path);
        size_t edge_length = g_.length(e);
        size_t path_length = CumulativeLength(g_, genome_path);
        DEBUG(
                "Diff length " << abs((int) edge_length - (int) path_length)
                        << "; genome path length " << path_length
                        << "; edge length " << edge_length);
        pair<string, pair<size_t, size_t>> c_id_and_pos = ContigIdAndPositions(
                e);

        if (c_id_and_pos.first == "")
            return;

        WriteAltPath(e, genome_path);
        size_t unique_violet = VioletLengthOfGenomeUnique(genome_path);
        map<TColorSet, size_t> color_cumm_lengths = ColorLengths(genome_path);
        if (color_cumm_lengths[kVioletColorSet] * 10
                > color_cumm_lengths[kBlueColorSet]) {
            WARN(
                    "Very long path along violet edges: "
                            << color_cumm_lengths[kVioletColorSet]
                            << " while blue path length: "
                            << color_cumm_lengths[kBlueColorSet]);
            WARN("While processing edge: " << g_.str(e));
        }
        if (color_cumm_lengths[kVioletColorSet] > 0)
            DEBUG("Violet edges in path");

        DEBUG("Total blue " << color_cumm_lengths[kBlueColorSet]);
        DEBUG("Total violet " << color_cumm_lengths[kVioletColorSet]);
        DEBUG("Total unique violet " << unique_violet);

        if (edge_length * path_length <= mem_lim) {
            size_t edit_dist = EditDistance(edge_nucls, path_nucls);
            DEBUG(
                    "Edit distance " << edit_dist << ". That is "
                            << double(edit_dist) / max(edge_length, path_length));
            pair<size_t, size_t> local_sim = LocalSimilarity(edge_nucls,
                    path_nucls);
            DEBUG(
                    "Local sim " << local_sim.first << " interval length "
                            << local_sim.second << " relative "
                            << ((double) local_sim.first / local_sim.second));
//            assembly_length-genome_length relative_local_sim genome_path_length assembly_length genome_length
//            contig_id contig_start contig_end genome_start genome_end min max local_sim sim_interval edit_dist edit_dist/max tot_blue
//            tot_violet unique_violet edge_id
            cerr
                    << str(
                            format(
                                    "%d %f %d %d %d %s %d %d %d %d %d %d %d %d %d %f %d %d %d %d")
                                    % ((int) edge_length - (int) path_length)
                                    % ((double) local_sim.first
                                            / local_sim.second)
                                    % genome_path.size() % edge_length
                                    % path_length % c_id_and_pos.first
                                    % c_id_and_pos.second.first
                                    % c_id_and_pos.second.second % genome_start
                                    % genome_end % min(edge_length, path_length)
                                    % max(edge_length, path_length)
                                    % local_sim.first % local_sim.second
                                    % edit_dist
                                    % (double(edit_dist)
                                            / max(edge_length, path_length))
                                    % color_cumm_lengths[kBlueColorSet]
                                    % color_cumm_lengths[kVioletColorSet]
                                    % unique_violet
                                    % g_.int_id(e)) << endl;
        } else {
            WARN("Edges were too long");
        }
    }

    void AnalyzeShortcutEdge(EdgeId e) {
        DEBUG("Analysing edge " << g_.str(e));
        pair<vector<EdgeId>, pair<size_t, size_t>> genome_path = FindGenomePath(
                g_.EdgeStart(e), g_.EdgeEnd(e), /*edge count bound*//*100*/300);
        if (!genome_path.first.empty()) {
            DEBUG(
                    "Non empty genome path of edge count "
                            << genome_path.first.size());
            DEBUG("Path " << g_.str(genome_path.first));
            Process(e, genome_path.first, genome_path.second.first,
                    genome_path.second.second);
        } else {
            DEBUG("Empty genome path");
        }
    }

public:
    SimpleInDelAnalyzer(const Graph& g, const ColorHandler<Graph>& coloring,
            const EdgesPositionHandler<Graph>& edge_pos,
            const vector<EdgeId> genome_path, TColorSet shortcut_color,
            const string& folder) :
            g_(g), coloring_(coloring), edge_pos_(edge_pos), genome_path_(
                    genome_path), shortcut_color_(shortcut_color), folder_(
                    folder) {
    }

    void Analyze() {
        cerr
                << "assembly_length-genome_length relative_local_sim genome_path_length assembly_length genome_length "
                << "contig_id contig_start contig_end genome_start genome_end min max local_sim sim_interval edit_dist "
                << "edit_dist/max tot_blue tot_violet unique_violet edge_id" << endl;
        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            if (coloring_.Color(*it) == shortcut_color_) {
                AnalyzeShortcutEdge(*it);
            }
        }
    }
private:
    DECL_LOGGER("SimpleInDelAnalyzer")
    ;
};

template<class gp_t>
class SimpleRearrangementDetector {
private:
    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const gp_t& gp_;
    const ColorHandler<Graph>& coloring_;
    const string ref_prefix_;
    const string folder_;
    mutable size_t count_;

    void ReportPossibleRearrangementConnection(EdgeId e, int start_ref_pos,
            int end_ref_pos, const string& folder) const {
        INFO(
                "Edge " << gp_.g.str(e)
                        << " identified as rearrangement connection");
        visualization::graph_labeler::LengthIdGraphLabeler<Graph> basic_labeler(gp_.g);
        visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(gp_.g, gp_.edge_pos);

        visualization::graph_labeler::CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);

        INFO(
                count_ << " example start_ref_pos: " << start_ref_pos
                        << " end_ref_pos: " << end_ref_pos);
        string filename = str(
                boost::format("%s%d_%d_%d_%d.dot") % folder % count_
                        % gp_.g.int_id(e) % start_ref_pos % end_ref_pos);
        GraphComponent<Graph> component = omnigraph::EdgeNeighborhood(gp_.g, e);
        visualization::visualization_utils::WriteComponent(component, filename, coloring_.ConstructColorer(component), labeler);
        count_++;
    }

    bool ContainsBlueEdge(const vector<EdgeId>& edges) const {
        for (size_t i = 0; i < edges.size(); ++i) {
            if (coloring_.Color(edges[i]) == kBlueColorSet)
                return true;
        }
        return false;
    }

    EdgeId GetBlueEdge(const vector<EdgeId>& edges) const {
        for (size_t i = 0; i < edges.size(); ++i) {
            if (coloring_.Color(edges[i]) == kBlueColorSet)
                return edges[i];
        }
        VERIFY(false);
        return EdgeId(NULL);
    }

    int GetRefPosition(EdgeId e, bool start_position) const {
        EdgePosition pos =
                RefPositions(gp_.edge_pos.GetEdgePositions(e)).front();
        int coeff = boost::ends_with(pos.contigId, "_RC") ? -1 : 1;
        Range range = pos.mr.initial_range;
        return coeff * (start_position ? range.start_pos : range.end_pos);
    }

    bool IsSingleRefPosition(EdgeId e) const {
        return RefPositions(gp_.edge_pos.GetEdgePositions(e)).size() == 1;
    }

    vector<EdgePosition> RefPositions(const vector<EdgePosition>& poss) const {
        vector < EdgePosition > answer;
        for (auto it = poss.begin(); it != poss.end(); ++it) {
            if (boost::starts_with(it->contigId, ref_prefix_)) {
                answer.push_back(*it);
            }
        }
        return answer;
    }

public:
    SimpleRearrangementDetector(const gp_t& gp,
            const ColorHandler<Graph>& coloring, const string& ref_prefix,
            const string& folder) :
            gp_(gp), coloring_(coloring), ref_prefix_(ref_prefix), folder_(
                    folder), count_(0) {
    }

    void Detect() const {
        for (auto it = gp_.g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            if (coloring_.Color(*it) == kRedColorSet) {
                INFO("Processing red edge " << gp_.g.str(*it));
                if (gp_.g.OutgoingEdgeCount(gp_.g.EdgeStart(*it)) == 2
                        && ContainsBlueEdge(
                                gp_.g.OutgoingEdges(gp_.g.EdgeStart(*it)))) {
                    EdgeId first_edge = GetBlueEdge(
                            gp_.g.OutgoingEdges(gp_.g.EdgeStart(*it)));
                    if (gp_.g.IncomingEdgeCount(gp_.g.EdgeEnd(*it)) == 2
                            && ContainsBlueEdge(
                                    gp_.g.IncomingEdges(gp_.g.EdgeEnd(*it)))) {
                        EdgeId second_edge = GetBlueEdge(
                                gp_.g.IncomingEdges(gp_.g.EdgeEnd(*it)));
                        if (first_edge != second_edge) {
                            INFO("Edges passed topology checks");
                            if (IsSingleRefPosition(first_edge)
                                    && IsSingleRefPosition(second_edge)) {
                                int start_ref_pos = GetRefPosition(first_edge,
                                        true);
                                int end_ref_pos = GetRefPosition(second_edge,
                                        false);
                                INFO("Edges had multiplicity one in reference");
                                ReportPossibleRearrangementConnection(*it,
                                        start_ref_pos, end_ref_pos, folder_);
                            } else {
                                INFO("Ooops");
                                INFO(
                                        "Edges had multiplicity more than one in reference");
                            }
                        }
                    }
                }
            }
        }
    }

    DECL_LOGGER("SimpleRearrangementDetector")
    ;
};

//template<class Graph>
//class GraphEdgeEnumerator {
//    const Graph& g_;
//    typedef typename Graph::EdgeId EdgeId;
//protected:
//    GraphEdgeEnumerator(const Graph& g) :
//            g_(g) {
//    }
//
//    const Graph& g() const {
//        return g_;
//    }
//public:
//    virtual ~GraphEdgeEnumerator() {
//    }
//    virtual map<EdgeId, string> Enumerate() const;
//};

//template<class Graph>
//class ThreadedGenomeEnumerator/*: public GraphEdgeEnumerator<Graph>*/{
////    typedef GraphEdgeEnumerator<Graph> base;
//    typedef typename Graph::EdgeId EdgeId;
//    const Graph& g_;
//    const vector<EdgeId> genome_path_;
//public:
//    ThreadedGenomeEnumerator(const Graph& g, const vector<EdgeId>& genome_path) :
//            g_(g), genome_path_(genome_path) {
//    }
//
//    /*virtual */
//    map<EdgeId, string> Enumerate() const {
//        map<EdgeId, string> answer;
//        //numerating genome path
//        int curr = 0;
//        for (auto it = genome_path_.begin(); it != genome_path_.end(); ++it) {
//            if (answer.find(*it) == answer.end()) {
//                curr++;
//                answer[*it] = ToString(curr);
//                answer[g_.conjugate(*it)] = ToString(-curr);
//            }
//        }
//        curr = 1000000;
//        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
//            if (answer.find(*it) == answer.end()) {
//                curr++;
//                answer[*it] = ToString(curr);
//                answer[g_.conjugate(*it)] = ToString(-curr);
//            }
//        }
//        return answer;
//    }
//};

//todo fixme use exact coordinates!
template<class Graph>
class BlockPrinter {

public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    BlockPrinter(const Graph& g, const CoordinatesHandler<Graph>& coords,
               const string& filename)
      : g_(g),
        coords_(coords),
        output_stream_(filename),
        curr_id_(1) {
      output_stream_
              << "genome_id\tcontig_name\tcanonical_id\tcontig_start_pos\tcontig_end_pos"
              << "\trefined_start_pos\trefined_end_pos\tsign\torig_id"
              << endl;
    }

    virtual ~BlockPrinter() {
    }

    //genome is supposed to perfectly correspond to some path in the graph
    void ProcessContig(unsigned genome_id, unsigned transparent_id,
                       const string& contig_name) {
        INFO("Processing contig " << transparent_id << " name " << contig_name);
        MappingPath<EdgeId> mapping_path = coords_.AsMappingPath(transparent_id);

        size_t graph_pos = 0;
        for (size_t i = 0; i < mapping_path.size(); ++i) {
            EdgeId e = mapping_path[i].first;
            MappingRange mapping = mapping_path[i].second;
            if (CheckPatternMatch(e)) {
                auto canon = CanonicalId(e);
                size_t next_graph_pos = graph_pos + g_.length(e);

                output_stream_
                        << (format("%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d")
                                % genome_id % contig_name % canon.first
                                % mapping.initial_range.start_pos
                                % mapping.initial_range.end_pos
                                % graph_pos
                                % next_graph_pos
                                % (canon.second ? "+" : "-") % g_.int_id(e)).str()
                        << endl;

                graph_pos = next_graph_pos;
            }
        }

//        VertexId v = g_.EdgeStart(coords_.FindGenomeFirstEdge(transparent_id));
//        size_t genome_pos = 0;
//
//        while (true) {
//            auto step = coords_.StepForward(v, transparent_id, genome_pos);
//            if (step.second == -1u)
//                break;
//
//            EdgeId e = step.first;
//
//            Range graph_pos(coords_.GetNewestPos(transparent_id, genome_pos),
//                            coords_.GetNewestPos(transparent_id, step.second));
//            Range contig_pos(
//                    coords_.GetOriginalPos(transparent_id, graph_pos.start_pos),
//                    coords_.GetOriginalPos(transparent_id, graph_pos.end_pos));
//            Range graph_pos_printable = coords_.GetPrintableRange(graph_pos);
//            Range contig_pos_printable = coords_.GetPrintableRange(contig_pos);
//
//            if (CheckPatternMatch(e)) {
//                auto canon = CanonicalId(e);
//
//                output_stream_
//                        << (format("%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d")
//                                % genome_id % contig_name % canon.first
//                                % contig_pos_printable.start_pos
//                                % contig_pos_printable.end_pos
//                                % graph_pos_printable.start_pos
//                                % graph_pos_printable.end_pos
//                                % (canon.second ? "+" : "-") % g_.int_id(e)).str()
//                        << endl;
//            }
//
//            v = g_.EdgeEnd(e);
//            genome_pos = step.second;
//        }
    }

    static void ConvertBlocksToGRIMM(const string &file_from,
                                     const string &file_to) {
        ifstream in(file_from);
        ofstream out(file_to);

        size_t id = 0;
        int last_genome_id = -1;
        size_t num_in_line = 0;
        while (!in.eof()) {
            ++id;

            string line;
            std::getline(in, line);
            if (id == 1)
                continue;

            if (line == "")
                continue;

            std::stringstream ss(line);

            int genome_id;
            string genome_name;
            string sign;
            size_t contig_id;

            string tmp;
            ss >> genome_id >> genome_name >> contig_id >> tmp >> tmp >> tmp
                    >> tmp >> sign >> tmp;
            if (genome_id != last_genome_id) {
                if (last_genome_id != -1)
                    out << "\n";
                out << "> " << genome_name << "\n";

                last_genome_id = genome_id;
                num_in_line = 0;
            }

            if (num_in_line > 10) {
                out << "\n";
                num_in_line = 0;
            }

            if (num_in_line != 0)
                out << " ";

            if (sign == "-")
                out << sign;
            out << contig_id;

            num_in_line++;
        }

        in.close();
        out.close();
    }

protected:
    virtual bool CheckPatternMatch(const EdgeId /* e */) {
        return true;
    }

    const Graph& g_;
    const CoordinatesHandler<Graph>& coords_;

private:
    ofstream output_stream_;
    size_t curr_id_;
    map<EdgeId, size_t> block_id_;

    pair<size_t, bool> CanonicalId(EdgeId e) {
//        size_t id = gp_.g.int_id(e);
//        size_t conj_id = gp_.g.int_id(gp_.g.conjugate(e));
//        if (id <= conj_id)
//            return make_pair(id, true);
//        else
//            return make_pair(conj_id, false);

        if (block_id_.count(e) > 0) {
            return make_pair(get(block_id_, e), true);
        } else if (block_id_.count(g_.conjugate(e)) > 0) {
            return make_pair(get(block_id_, g_.conjugate(e)), false);
        } else {
            block_id_[e] = curr_id_++;
            return make_pair(get(block_id_, e), true);
        }
    }

    DECL_LOGGER("BlockPrinter");
};

template<class Graph>
class UniqueBlockPrinter : public BlockPrinter<Graph> {
 public:
  UniqueBlockPrinter(const Graph &g, const CoordinatesHandler<Graph> &coords,
            const string &filename, const vector<pair<size_t, size_t>> rc_pairs)
      : BlockPrinter<Graph>(g, coords, filename),
        rc_pairs_(rc_pairs),
        contig_map_(),
        cur_time_(rc_pairs_.size()),
        glob_time_(0) {
    PrepareContigMap(rc_pairs);
  }

  // virtual ~UniqueBlockPrinter() {
  // }

 protected:
  virtual bool CheckPatternMatch(const EdgeId e) {
    glob_time_++;

    const auto &ranges = this->coords_.GetRawRanges(e);
    if (ranges.size() != rc_pairs_.size())
      return false;

    for (const auto &e : ranges) {
      size_t my_id = contig_map_.at(e.first);
      if (cur_time_[my_id] == glob_time_)
        return false;
      cur_time_[my_id] = glob_time_;
    }

    // By the Dirichlet priciple..
    return true;
  }

 private:
  void PrepareContigMap(const vector<pair<size_t, size_t>> rc_pairs) {
    for (size_t i = 0; i < rc_pairs.size(); ++i) {
      const auto &p = rc_pairs[i];

      contig_map_[p.first] = i;
      contig_map_[p.second] = i;
    }
  }

  vector<pair<size_t, size_t>> rc_pairs_;
  unordered_map<size_t, size_t> contig_map_;
  vector<size_t> cur_time_;
  size_t glob_time_;

  DECL_LOGGER("UniqueBlockPrinter");
};

//template<class Graph, class Mapper>
//class ContigBlockStats {
//    typedef ThreadedGenomeEnumerator<Graph> Enumerator;
//    const Graph& g_;
//    const EdgesPositionHandler<Graph>& edge_pos_;
//    const Mapper mapper_;
//    const vector<EdgeId> genome_path_;
//    ContigStream& contigs_;
//    const map<EdgeId, string> labels_;
//
//    const string& get(const map<EdgeId, string>& from, EdgeId key) const {
//        auto it = from.find(key);
//        VERIFY(it != from.end());
//        return it->second;
//    }
//
//    void ReportGenomeBlocks() const {
//        set < EdgeId > visited;
////        cerr << "Genome blocks started" << endl;
//        for (auto it = genome_path_.begin(); it != genome_path_.end(); ++it) {
//            if (visited.count(*it) > 0)
//                continue;
//            cerr << get(labels_, *it) << " $ "
//                    << g_.int_id(*it)
//                    << " $ "
//                    << g_.length(*it)
//                    /*<< " positions: " << edge_pos_.GetEdgePositions(*it) */<< endl;
//            visited.insert(*it);
//            visited.insert(g_.conjugate(*it));
//        }
////        cerr << "Genome blocks ended" << endl;
//    }
//
//    void ReportOtherBlocks() const {
////        cerr << "Other blocks started" << endl;
//        for (auto it = labels_.begin(); it != labels_.end(); ++it) {
//            if (boost::lexical_cast<int>(it->second) > (int) 1000000) {
//                cerr << get(labels_, it->first) << " $ "
//                        << g_.int_id(it->first)
//                        << " $ "
//                        << g_.length(it->first)
//                        /*<< " positions: " << edge_pos_.GetEdgePositions(it->first) */<< endl;
//            }
//        }
////        cerr << "Other blocks ended" << endl;
//    }
//
//    void ReportContigs() const {
//        contigs_.reset();
//        Contig contig;
//        cerr << "Contigs started" << endl;
//        while (!contigs_.eof()) {
//            contigs_ >> contig;
//            vector < EdgeId > path =
//                    mapper_.MapSequence(contig.sequence()).simple_path();
//            cerr << contig.name() << " $ ";
//            string delim = "";
//            for (auto it = path.begin(); it != path.end(); ++it) {
//                cerr << delim << get(labels_, *it);
//                delim = ";";
//            }
//            cerr << endl;
//        }
//        cerr << "Contigs ended" << endl;
//    }
//
//public:
//    ContigBlockStats(const Graph& g,
//            const EdgesPositionHandler<Graph>& edge_pos, const Mapper& mapper,
//            const Sequence& genome, ContigStream& contigs) :
//            g_(g), edge_pos_(edge_pos), mapper_(mapper), genome_path_(
//                    mapper_.MapSequence(genome).simple_path()), contigs_(
//                    contigs), labels_(Enumerator(g_, genome_path_).Enumerate()) {
//    }
//
//    void Count() const {
//        cerr << "Block id $ Graph edge id $ (for debug) $ Length (in 201-mers)"
//                << endl;
//
//        ReportGenomeBlocks();
//        ReportOtherBlocks();
//
//        cerr << "Contig id $ Block ids" << endl;
//        ReportContigs();
//    }
//};

template<class Graph>
class AlternatingPathsCounter {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph& g_;
    const ColorHandler<Graph>& coloring_;

    TColorSet InvertColor(TColorSet color) const {
        if (color == kRedColorSet) {
            return kBlueColorSet;
        } else if (color == kBlueColorSet) {
            return kRedColorSet;
        }
        VERIFY(false);
        return kBlueColorSet;
    }

    vector<EdgeId> FilterEdges(vector<EdgeId> edges, TColorSet color) const {
        vector<EdgeId> answer;
        for (size_t i = 0; i < edges.size(); ++i) {
            if (coloring_.Color(edges[i]) == color) {
                answer.push_back(edges[i]);
            }
        }
        return answer;
    }

    vector<EdgeId> OutgoingEdges(VertexId v, TColorSet color) const {
        DEBUG(
                "Looking for outgoing edges for vertex " << g_.str(v)
                        << " of color " << color);
        return FilterEdges(g_.OutgoingEdges(v), color);
    }

    vector<EdgeId> IncomingEdges(VertexId v, TColorSet color) const {
        DEBUG(
                "Looking for incoming edges for vertex " << g_.str(v)
                        << " of color " << color);
        return FilterEdges(g_.IncomingEdges(v), color);
    }

    bool CheckNotContains(vector<EdgeId>& path, EdgeId e) const {
        return std::find(path.begin(), path.end(), e) == path.end();
    }

    VertexId OtherVertex(EdgeId e, VertexId v) const {
        VERIFY(
                g_.EdgeStart(e) != g_.EdgeEnd(e)
                        && (g_.EdgeStart(e) == v || g_.EdgeEnd(e) == v));
        if (g_.EdgeStart(e) == v) {
            DEBUG("Next vertex " << g_.EdgeEnd(e));
            return g_.EdgeEnd(e);
        }
        DEBUG("Next vertex " << g_.EdgeStart(e));
        return g_.EdgeStart(e);
    }

    bool Grow(vector<EdgeId>& path, VertexId last_vertex) const {
        DEBUG("Growing path for vertex " << g_.str(last_vertex));
        EdgeId last_edge = path.back();
        DEBUG("Last edge " << last_edge);
        TColorSet next_color = InvertColor(coloring_.Color(last_edge));
        vector<EdgeId> next_candidates =
                (g_.EdgeEnd(last_edge) == last_vertex) ?
                        IncomingEdges(last_vertex, next_color) :
                        OutgoingEdges(last_vertex, next_color);
        if (next_candidates.empty()) {
            DEBUG("No candidates");
            return true;
        }
        if (next_candidates.size() > 1) {
            DEBUG("Several candidates");
            return false;
        }
        EdgeId next_edge = next_candidates.front();
        DEBUG(
                "Adding edge " << g_.str(next_edge) << " of color "
                        << coloring_.Color(next_edge));
        if (!CheckNotContains(path, next_edge)) {
            WARN("PROBLEM");
            return false;
        }

        path.push_back(next_edge);
        return Grow(path, OtherVertex(next_edge, last_vertex));
    }

    vector<EdgeId> AlternatingPathContainingEdge(EdgeId e) const {
        vector<EdgeId> answer;
        vector<EdgeId> tmp_path;
        DEBUG("Growing backward");
        tmp_path.push_back(e);
        if (Grow(tmp_path, g_.EdgeStart(e))) {
            answer.insert(answer.end(), tmp_path.rbegin(), tmp_path.rend());
            tmp_path.clear();
            DEBUG("Growing forward");
            tmp_path.push_back(e);
            if (Grow(tmp_path, g_.EdgeEnd(e))) {
                answer.insert(answer.end(), (++tmp_path.begin()),
                        tmp_path.end());
                return answer;
            }
        }
        return vector<EdgeId>();
    }

    void ProcessAltPath(const vector<EdgeId>& path) const {
        DEBUG("Processing path of length " << path.size());
        cerr << path.size() << endl;
    }

public:
    AlternatingPathsCounter(const Graph& g, const ColorHandler<Graph>& coloring) :
            g_(g), coloring_(coloring) {
    }

    void CountPaths() const {
        set<EdgeId> visited_edges;
        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            if (visited_edges.count(*it) > 0)
                continue;
            if (coloring_.Color(*it) == kRedColorSet) {
                DEBUG("Looking for alt path for edge " << g_.str(*it));
                vector<EdgeId> alt_path = AlternatingPathContainingEdge(*it);
                if (!alt_path.empty()) {
                    ProcessAltPath(alt_path);
                    visited_edges.insert(alt_path.begin(), alt_path.end());
                }
            }
        }
    }
private:
    DECL_LOGGER("AlternatingPathsCounter")
    ;
};

template<class Graph, class Mapper>
class MissingGenesAnalyser {
    typedef typename Graph::EdgeId EdgeId;
    const Graph& g_;
    const ColorHandler<Graph>& coloring_;
    const EdgesPositionHandler<Graph>& edge_pos_;
    const Sequence genome_;
    const Mapper mapper_;
    const vector<pair<bool, pair<size_t, size_t>>> locations_;
    const string output_dir_;

    void ReportLocality(const Sequence& s, const string& out_file) {
        visualization::graph_labeler::LengthIdGraphLabeler<Graph> basic_labeler(g_);
        visualization::graph_labeler::EdgePosGraphLabeler<Graph> pos_labeler(g_, edge_pos_);

        visualization::graph_labeler::CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);

        WriteComponentsAlongPath(g_, labeler, out_file, /*split_length*/1000, /*vertex_number*/15
                , mapper_.MapSequence(s), *ConstructBorderColorer(g_, coloring_));
    }

public:
    MissingGenesAnalyser(const Graph& g, const ColorHandler<Graph>& coloring
            , const EdgesPositionHandler<Graph>& edge_pos
            , const Sequence& genome
            , const Mapper& mapper
            , const vector<pair<bool, pair<size_t, size_t>>>& locations
            , const string& output_dir):
    g_(g), coloring_(coloring), edge_pos_(edge_pos), genome_(genome), mapper_(mapper), locations_(locations), output_dir_(output_dir) {

    }

    void Analyze() {
        remove_dir(output_dir_);
        make_dir(output_dir_);
        for (size_t i = 0; i < locations_.size(); ++i) {
            pair<bool, pair<size_t, size_t>> location = locations_[i];
            Sequence locality = genome_.Subseq(location.second.first - g_.k(), location.second.second + g_.k());
            if (location.first) {
                locality = !locality;
            }
            ReportLocality(locality, output_dir_ + ToString(i) + ".dot");
        }
    }
};}
