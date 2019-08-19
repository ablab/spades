#pragma once

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/components/graph_component.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "visualization/printing_parameter_storage.hpp"
//#include "edges_position_handler.hpp"

using namespace omnigraph;

namespace visualization {

namespace graph_colorer {

template<typename ElementId>
class ElementColorer : public virtual printing_parameter_storage::ParameterStorage<ElementId, std::string> {
public:
    template<typename Iter>
    std::set<ElementId> ColoredWith(Iter begin, Iter end, const std::string &color) {
        std::set<ElementId> result;
        for (Iter it = begin; it != end; ++it) {
            if (this->GetValue(*it) == color)
                result.insert(*it);
        }
        return result;
    }
};

//TODO remove all default color parameters!

template<typename ElementId>
class MapColorer : public ElementColorer<ElementId>, public printing_parameter_storage::MapParameterStorage<ElementId, std::string> {
public:
    MapColorer(const std::string &default_color) : printing_parameter_storage::MapParameterStorage<ElementId, std::string>(default_color) {
    }

    MapColorer(const std::map<ElementId, std::string> &color_map) : printing_parameter_storage::MapParameterStorage<ElementId, std::string>(color_map) {
    }

    MapColorer(const std::map<ElementId, std::string> &color_map, const std::string &default_color)
            : printing_parameter_storage::MapParameterStorage<ElementId, std::string>(color_map, default_color) {
    }

    template<class It>
    MapColorer(It begin, It end, const std::string &color, const std::string &default_color)
            : printing_parameter_storage::MapParameterStorage<ElementId, std::string>(begin, end, color, default_color) {
    }

    virtual ~MapColorer() {
    }
};

template<typename ElementId>
class FixedColorer : public MapColorer<ElementId> {
public:
    FixedColorer(const std::string &default_color) : MapColorer<ElementId>(default_color) {
    }
};

template<class Graph>
class SetColorer : public MapColorer<typename Graph::EdgeId> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;

    template<class It>
    std::map<EdgeId, std::string> ConstructColorMap(It begin, It end, const std::string &color) {
        std::map<EdgeId, std::string> result;
        for (auto it = begin; it != end; ++it) {
            result[*it] = color;
        }
        return result;
    }

public:
    template<class It>
    SetColorer(const Graph &graph, It begin, It end, const std::string &color) :
            MapColorer<typename Graph::EdgeId>(ConstructColorMap(begin, end, color), "black"), graph_(graph) {
    }

    template<class Collection>
    SetColorer(const Graph &graph, const Collection &c, const std::string &color) :
            MapColorer<typename Graph::EdgeId>(ConstructColorMap(c.begin(), c.end(), color), "black"),
            graph_(graph) {
    }

};
//
//template<class Graph>
//class PositionsEdgeColorer: public ElementColorer<typename Graph::EdgeId> {
//private:
//    typedef typename Graph::VertexId VertexId;
//    typedef typename Graph::EdgeId EdgeId;
//    const Graph &graph_;
//    EdgesPositionHandler<Graph> &positions_;
//public:
//    PositionsEdgeColorer(const Graph &graph, EdgesPositionHandler<Graph> &positions):
//            graph_(graph), positions_(positions)  {
//    }
//    string GetValue(EdgeId element) const {
//        std::vector<EdgeId> path;
//        path.push_back(element);
//        if (positions_.GetEdgePositions(element).size() == 0) return "black";
//        else {
//            if (positions_.IsConsistentWithGenome(path)) return "green";
//            else return "orange";
//        }
//    }
//};


template<class Graph>
class CompositeEdgeColorer : public ElementColorer<typename Graph::EdgeId> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    std::string default_color_;
    std::vector<std::shared_ptr<ElementColorer<typename Graph::EdgeId>>> colorers_;

    std::vector<std::string> CollectColors(EdgeId edge) const {
        std::vector<std::string> result = {default_color_};
        for (auto it = colorers_.begin(); it != colorers_.end(); ++it) {
            std::string next_color = (*it)->GetValue(edge);
            if (std::find(result.begin(), result.end(), next_color) == result.end())
                result.push_back(next_color);
        }
        return result;
    }

    std::string ConstructColorString(const std::vector<std::string> &colors) const {
        if (colors.size() == 1)
            return default_color_;
        std::string result;
        for (size_t i = 1; i < colors.size(); i++)
            result += ":" + colors[i];
        return result.substr(1, result.size());
    }

public:
    CompositeEdgeColorer(const std::string &default_color) : default_color_(default_color) {
    }

    CompositeEdgeColorer(std::shared_ptr<ElementColorer<typename Graph::EdgeId>> colorer,
                         const std::string &default_color) : default_color_(default_color) {
        AddColorer(colorer);
    }

    CompositeEdgeColorer(std::shared_ptr<ElementColorer<typename Graph::EdgeId>> colorer1,
                         std::shared_ptr<ElementColorer<typename Graph::EdgeId>> colorer2,
                         const std::string &default_color) : default_color_(default_color) {
        AddColorer(colorer1);
        AddColorer(colorer2);
    }

    void AddColorer(std::shared_ptr<ElementColorer<typename Graph::EdgeId>> colorer) {
        colorers_.push_back(colorer);
    }

    void PopColorer() {
        colorers_.pop_back();
    }

    std::string GetValue(EdgeId edge) const {
        return ConstructColorString(CollectColors(edge));
    }
};

template<class Graph>
class GraphColorer
        : public ElementColorer<typename Graph::VertexId>, public ElementColorer<typename Graph::EdgeId> {
public:
    std::string GetValue(typename Graph::VertexId) const = 0;

    std::string GetValue(typename Graph::EdgeId) const = 0;

    template<typename Iter>
    std::set<typename Iter::value_type> ColoredWith(Iter begin, Iter end, const std::string &color) {
        return ElementColorer<typename Iter::value_type>::ColoredWith(begin, end, color);
    }
};

template<class Graph>
class DelegatingGraphColorer : public GraphColorer<Graph> {
private:
    const GraphColorer<Graph> &inner_colorer_;
public:
    DelegatingGraphColorer(const GraphColorer<Graph> &inner_colorer) : inner_colorer_(inner_colorer) {
    }

    std::string GetValue(typename Graph::VertexId v) const {
        return inner_colorer_.GetValue(v);
    }

    std::string GetValue(typename Graph::EdgeId e) const {
        return inner_colorer_.GetValue(e);
    }
};

template<typename Graph>
class BorderDecorator : public GraphColorer<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const GraphComponent<Graph> &component_;
//    const shared_ptr<const ElementColorer<typename Graph::VertexId>> vertex_colorer_ptr_;
//    const shared_ptr<const ElementColorer<typename Graph::EdgeId>> edge_colorer_ptr_;
    const ElementColorer<typename Graph::VertexId> &vertex_colorer_;
    const ElementColorer<typename Graph::EdgeId> &edge_colorer_;
    const std::string border_color_;
public:
//    BorderDecorator(const GraphComponent<Graph> &component,
//            const shared_ptr<const GraphColorer<Graph>> colorer,
//            const string &border_color) :
//            component_(component), vertex_colorer_ptr_(colorer), edge_colorer_ptr_(
//                    colorer), vertex_colorer_(*colorer), edge_colorer_(
//                    *colorer), border_color_(border_color) {
//    }

    BorderDecorator(const GraphComponent<Graph> &component, const GraphColorer<Graph> &colorer,
                    const std::string &border_color = "yellow") :
            component_(component), vertex_colorer_(colorer), edge_colorer_(colorer),
            border_color_(border_color) {
    }

    std::string GetValue(VertexId v) const {
        if (component_.IsBorder(v)) {
            return border_color_;
        } else {
            return vertex_colorer_.GetValue(v);
        }
    }

    std::string GetValue(EdgeId e) const {
        return edge_colorer_.GetValue(e);
    }

    static std::shared_ptr<BorderDecorator<Graph>> GetInstance(const GraphComponent<Graph> &component,
                                                               const GraphColorer<Graph> &colorer,
                                                               const std::string &border_color = "yellow") {
        return std::make_shared<BorderDecorator<Graph>>(component, colorer, border_color);
    }
};


template<typename Graph>
class SinkSourceDecorator : public GraphColorer<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const GraphComponent<Graph> &component_;
//  const shared_ptr<const ElementColorer<typename Graph::VertexId>> vertex_colorer_ptr_;
//  const shared_ptr<const ElementColorer<typename Graph::EdgeId>> edge_colorer_ptr_;
    const ElementColorer<typename Graph::VertexId> &vertex_colorer_;
    const ElementColorer<typename Graph::EdgeId> &edge_colorer_;
    const std::string sink_color_, source_color_, sinksource_color_;
public:

    SinkSourceDecorator(const GraphComponent<Graph> &component,
                        const GraphColorer<Graph> &colorer, const std::string &sink_color = "red",
                        const std::string &source_color = "orange", const std::string &sinksource_color = "green") :
            component_(component), vertex_colorer_(colorer), edge_colorer_(colorer), sink_color_(sink_color),
            source_color_(source_color), sinksource_color_(sinksource_color) {
    }

    std::string GetValue(VertexId v) const {
        if (component_.exits().count(v) && !component_.entrances().count(v)) {
            return sink_color_;
        }
        if (component_.entrances().count(v) && !component_.exits().count(v)) {
            return source_color_;
        }
        if (component_.entrances().count(v) && component_.exits().count(v)) {
            return sinksource_color_;
        }

        return vertex_colorer_.GetValue(v);
    }

    std::string GetValue(EdgeId e) const {
        return edge_colorer_.GetValue(e);
    }

    static std::shared_ptr<SinkSourceDecorator<Graph>> GetInstance(const GraphComponent<Graph> &component,
                                                                   const GraphColorer<Graph> &colorer,
                                                                   const std::string &sink_color = "red",
                                                                   const std::string &source_color = "orange") {
        return std::make_shared<SinkSourceDecorator<Graph>>(component, colorer, sink_color, source_color);
    }
};

template<class Graph>
class CompositeGraphColorer : public GraphColorer<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const std::shared_ptr<ElementColorer<VertexId>> vertex_colorer_;
    const std::shared_ptr<ElementColorer<EdgeId>> edge_colorer_;
public:
    CompositeGraphColorer(std::shared_ptr<ElementColorer<VertexId>> vertex_colorer,
                          std::shared_ptr<ElementColorer<EdgeId>> edge_colorer) :
            vertex_colorer_(vertex_colorer),
            edge_colorer_(edge_colorer) {
    }

//    explicit CompositeGraphColorer(shared_ptr<ElementColorer<EdgeId>> edge_colorer = make_shared<FixedColorer<EdgeId>>("black")) :
//                vertex_colorer_(shared_ptr<ElementColorer<VertexId>>(new FixedColorer<VertexId>("white"))),
//                edge_colorer_(edge_colorer) {
//    }

    std::string GetValue(VertexId v) const {
        return vertex_colorer_->GetValue(v);
    }

    std::string GetValue(EdgeId e) const {
        return edge_colorer_->GetValue(e);
    }

};


// edge_colorer management is passed here
//TODO check all usages
template<class Graph>
std::shared_ptr<GraphColorer<Graph>> DefaultColorer(const Graph & /*g*/,
                                                    std::shared_ptr<ElementColorer<typename Graph::EdgeId>> edge_colorer) {
    return std::shared_ptr<GraphColorer<Graph>>(
            new CompositeGraphColorer<Graph>(std::make_shared<FixedColorer<typename Graph::VertexId>>("white"),
                                             edge_colorer));
}

template<class Graph>
std::shared_ptr<GraphColorer<Graph>> DefaultColorer(const Graph &g,
                                                    const Path<typename Graph::EdgeId> &path1,
                                                    const Path<typename Graph::EdgeId> &path2) {
    std::shared_ptr<ElementColorer<typename Graph::EdgeId>> edge_colorer =
            std::make_shared<CompositeEdgeColorer<Graph>>(
                    std::make_shared<SetColorer<Graph>>(g, path1.sequence(), "red"),
                    std::make_shared<SetColorer<Graph>>(g, path2.sequence(), "blue"), "black");
    return DefaultColorer(g, edge_colorer);
}

template<class Graph>
std::shared_ptr<GraphColorer<Graph>> DefaultColorer(const Graph & /*g*/) {
    return std::shared_ptr<GraphColorer<Graph>>(new CompositeGraphColorer<Graph>(
            std::make_shared<FixedColorer<typename Graph::VertexId>>("white"),
            std::make_shared<FixedColorer<typename Graph::EdgeId>>("black")));
}
}
}
