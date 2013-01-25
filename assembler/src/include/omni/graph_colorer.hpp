//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "omni_utils.hpp"
#include "edges_position_handler.hpp"

namespace omnigraph {

template<typename ElementId>
class ElementColorer {
public:
	virtual string GetColour(ElementId element) const = 0;

	template<class It>
	map<ElementId, string> GetColours(It begin, It end) const {
		return GetColours(set<ElementId>(begin, end));
	}

	virtual map<ElementId, string> GetColours(const set<ElementId> &elements) const {
		map<ElementId, string> result;
		for(auto it = elements.begin(); it != elements.end(); ++it) {
			result[*it] = GetColour(*it);
		}
		return result;
	}

	virtual ~ElementColorer() {
	}
};

template<typename ElementId>
class FixedColorer: public ElementColorer<ElementId> {
	string default_color_;
public:
	FixedColorer(const string& default_color):
		default_color_(default_color) {

	}

	string GetColour(ElementId element) const {
		return default_color_;
	}

	virtual ~FixedColorer() {
	}
};



//template<class Graph>
//class VertexColorer : public ElementColorer<Graph, typename Graph::VertexId> {
//
//};
//
//template<class Graph>
//class EdgeColorer : public ElementColorer<Graph, typename Graph::EdgeId> {
//
//};
//
template<class Graph>
class GraphColorer {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:

	virtual string GetColour(VertexId v) const = 0;

	virtual map<VertexId, string> GetColours(const set<VertexId> &vertices) const = 0;

	virtual string GetColour(EdgeId e) const = 0;

	virtual map<EdgeId, string> GetColours(const set<EdgeId> &edges) const = 0;

	virtual ~GraphColorer() {
	}
};

template<class Graph>
class CompositeGraphColorer: public GraphColorer<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	const auto_ptr<ElementColorer<VertexId>> vertex_colorer_;
	const auto_ptr<ElementColorer<EdgeId>> edge_colorer_;
public:
//	CompositeGraphColorer(auto_ptr<ElementColorer<Graph, VertexId>> vertex_colorer
//			, auto_ptr<ElementColorer<Graph, EdgeId>> edge_colorer) :
//				vertex_colorer_(vertex_colorer),
//				edge_colorer_(edge_colorer) {
//
//	}

	CompositeGraphColorer(ElementColorer<VertexId>* vertex_colorer
			, ElementColorer<EdgeId>* edge_colorer) :
				vertex_colorer_(vertex_colorer),
				edge_colorer_(edge_colorer) {

	}

	/*virtual */string GetColour(VertexId v) const {
		return vertex_colorer_->GetColour(v);
	}

	/*virtual */map<VertexId, string> GetColours(const set<VertexId> &vertices) const {
		return vertex_colorer_->GetColours(vertices.begin(), vertices.end());
	}

	/*virtual */string GetColour(EdgeId e) const {
		return edge_colorer_->GetColour(e);
	}

	/*virtual */map<EdgeId, string> GetColours(const set<EdgeId> &edges) const {
		return edge_colorer_->GetColours(edges.begin(), edges.end());
	}

};

template<typename ElementId>
class MapColorer : public ElementColorer<ElementId> {
private:
	map<ElementId, string> color_map_;
	optional<string> default_color_;
public:
	MapColorer(const map<ElementId, string> &color_map) : color_map_(color_map) {
	}

	MapColorer(const map<ElementId, string> &color_map, const string& default_color) :
		color_map_(color_map),
		default_color_(default_color) {
	}

	template<class It>
	MapColorer(It begin, It end, const string& color, const string& default_color) :
		default_color_(default_color) {
		for (auto it = begin; it != end; ++it) {
			color_map_.insert(make_pair(*it, color));
		}
	}

	virtual ~MapColorer() {
	}

	virtual string GetColour(ElementId element) const {
		if (color_map_.count(element) != 0) {
			return color_map_.find(element)->second;
//			return color_map_[element];
		} else {
			if (default_color_) {
				return *default_color_;
			} else {
				VERIFY(false);
				return "";
			}
		}
	}

};

template<class Graph>
class BorderVertexColorer : public ElementColorer<typename Graph::VertexId> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;

	bool IsBorder(VertexId v, const set<VertexId> &vs) const {
		const vector<EdgeId> outgoing_edges = graph_.OutgoingEdges(v);
		const vector<EdgeId> incoming_edges = graph_.IncomingEdges(v);
		set<EdgeId> adjacent_edges;
		adjacent_edges.insert(outgoing_edges.begin(), outgoing_edges.end());
		adjacent_edges.insert(incoming_edges.begin(), incoming_edges.end());
		for (auto e_it = adjacent_edges.begin(); e_it != adjacent_edges.end();
				++e_it) {
			if (vs.count(graph_.EdgeStart(*e_it)) == 0
					|| vs.count(graph_.EdgeEnd(*e_it)) == 0) {
				return true;
			}
		}
		return false;
	}

public:

	BorderVertexColorer(const Graph &graph) :
			graph_(graph) {
	}

	virtual string GetColour(VertexId element) const {
		//TODO Explain why
		//VERIFY(false);
		return "white";
	}

	virtual map<VertexId, string> GetColours(const set<VertexId> &elements) const {
		map<VertexId, string> result;
		for(auto it = elements.begin(); it != elements.end(); ++it) {
			string value = IsBorder(*it, elements) ? "yellow" : "white";
			result[*it] = value;
		}
		return result;
	}
};

template<class Graph>
class PathColorer {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	const Path<EdgeId> &path1_;
	const Path<EdgeId> &path2_;

	void SetColor(map<EdgeId, string> &color, EdgeId edge, string col) const {
		if (color.count(edge) != 0 && color[edge] != col) {
			color[edge] = "purple";
		} else
			color[edge] = col;
	}

	void ConstructColorMap(map<EdgeId, string> &color) const {
		for (auto it = path1_.sequence().begin(); it != path1_.sequence().end();
				++it) {
			SetColor(color, *it, "red");
		}
		for (auto it = path2_.sequence().begin(); it != path2_.sequence().end();
				++it) {
			SetColor(color, *it, "blue");
		}
	}

	void ConstructBlackEdgesSet(set<EdgeId> &result) const {
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd();
				++iterator) {
			result.insert(*iterator);
		}
		for (auto iterator = path1_.sequence().begin();
				iterator != path1_.sequence().end(); ++iterator) {
			result.erase(*iterator);
		}
		for (auto iterator = path2_.sequence().begin();
				iterator != path2_.sequence().end(); ++iterator) {
			result.erase(*iterator);
		}
	}

public:
	PathColorer(const Graph &graph, const Path<EdgeId> &path1,
			const Path<EdgeId> &path2) :
			graph_(graph), path1_(path1), path2_(path2) {
	}

	map<EdgeId, string> ColorPath() const {
		map<EdgeId, string> colors;
		ConstructColorMap(colors);
		return colors;
	}

	set<EdgeId> BlackEdges() const {
		set<EdgeId> result;
		ConstructBlackEdgesSet(result);
		return result;
	}
};

template<class Graph>
class NewPathColorer {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;

	void SetColor(map<EdgeId, string> &color, EdgeId edge, string col) const {
		if (color.count(edge) != 0 && color[edge] != col) {
			color[edge] = "purple";
		} else
			color[edge] = col;
	}

	void ConstructColorMap(map<EdgeId, string> &color, const Path<EdgeId> &path1,
			const Path<EdgeId> &path2) const {
		for (auto it = path1.sequence().begin(); it != path1.sequence().end();
				++it) {
			SetColor(color, *it, "red");
		}
		for (auto it = path2.sequence().begin(); it != path2.sequence().end();
				++it) {
			SetColor(color, *it, "blue");
		}
	}

	void ConstructBlackEdgesSet(set<EdgeId> &result, const Path<EdgeId> &path1,
			const Path<EdgeId> &path2) const {
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd();
				++iterator) {
			result.insert(*iterator);
		}
		for (auto iterator = path1.sequence().begin();
				iterator != path1.sequence().end(); ++iterator) {
			result.erase(*iterator);
		}
		for (auto iterator = path2.sequence().begin();
				iterator != path2.sequence().end(); ++iterator) {
			result.erase(*iterator);
		}
	}

public:
	NewPathColorer(const Graph &graph) :
			graph_(graph) {
	}

	map<EdgeId, string> ColorPath(const Path<EdgeId> &path1,
			const Path<EdgeId> &path2) const {
		map<EdgeId, string> colors;
		ConstructColorMap(colors, path1, path2);
		return colors;
	}

	set<EdgeId> BlackEdges(const Path<EdgeId> &path1,
			const Path<EdgeId> &path2) const {
		set<EdgeId> result;
		ConstructBlackEdgesSet(result, path1, path2);
		return result;
	}
};

// edge_colorer management is passed here
template <class Graph>
auto_ptr<GraphColorer<Graph>> DefaultColorer(const Graph& g,
		ElementColorer<typename Graph::EdgeId>* edge_colorer) {
	return auto_ptr<GraphColorer<Graph>>(new CompositeGraphColorer<Graph>(new BorderVertexColorer<Graph>(g), edge_colorer));
}

template <class Graph>
auto_ptr<GraphColorer<Graph>> DefaultColorer(const Graph& g,
		const map<typename Graph::EdgeId, string>& edge_color_map,
		const string& default_color = "") {
	return DefaultColorer(g, new MapColorer<typename Graph::EdgeId>(edge_color_map, default_color));
}

template <class Graph>
auto_ptr<GraphColorer<Graph>> DefaultColorer(const Graph& g,
		const Path<typename Graph::EdgeId>& path1,
		const Path<typename Graph::EdgeId>& path2) {
	return DefaultColorer(g, PathColorer<Graph>(g, path1, path2).ColorPath());
}

template <class Graph>
auto_ptr<GraphColorer<Graph>> DefaultColorer(const Graph& g) {
	map<typename Graph::EdgeId, string> empty_map;
	return DefaultColorer(g, empty_map);
}






template<class Graph>
class PositionsEdgeColorer: public ElementColorer<typename Graph::EdgeId> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	EdgesPositionHandler<Graph> &positions_;
public:
	PositionsEdgeColorer(const Graph &graph, EdgesPositionHandler<Graph> &positions):
			graph_(graph), positions_(positions)  {
	}
	virtual string GetColour(EdgeId element) const {
		std::vector<EdgeId> path;
		path.push_back(element);
		if (positions_.GetEdgePositions(element).size() == 0) return "black";
		else {
			if (positions_.IsConsistentWithGenome(path)) return "green";
			else return "orange";
		}
	}

};

}
