#pragma once

namespace compare {

enum edge_type {
	//don't change order!!!
	black = 0,
	red,
	blue,
	violet
};

template<class Graph, class Element>
class ElementColorHandler: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;

	restricted::map<Element, edge_type> data_;
public:
	static string color_str(edge_type color) {
		static string colors[] = { "black", "red", "blue", "purple" };
		return colors[(int) color];
	}

	ElementColorHandler(const Graph& g) :
			base(g, "ElementColorHandler") {

	}

	void Paint(Element e, edge_type color) {
		data_[e] = (edge_type) ((int) data_[e] | (int) color);
	}

	edge_type Color(Element e) const {
		auto it = data_.find(e);
		if (it == data_.end())
			return edge_type::black;
		else
			return it->second;
	}

	string ColorStr(Element e) const {
		return color_str(Color(e));
	}

	/*virtual*/
	void HandleDelete(Element e) {
		data_.erase(e);
	}

};

template<class Graph>
class ColorHandler: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	ElementColorHandler<Graph, EdgeId> edge_color_;
	ElementColorHandler<Graph, VertexId> vertex_color_;
public:

	ColorHandler(const Graph& g) :
			base(g, "ColorHandler"), edge_color_(g), vertex_color_(g) {

	}

	void Paint(EdgeId e, edge_type color) {
		edge_color_.Paint(e, color);
	}

	void Paint(VertexId v, edge_type color) {
		vertex_color_.Paint(v, color);
	}

	edge_type Color(EdgeId e) const {
		return edge_color_.Color(e);
	}

	edge_type Color(VertexId v) const {
		return vertex_color_.Color(v);
	}

	map<EdgeId, string> EdgeColorMap() const {
		map<EdgeId, string> answer;
		for (auto it = this->g().SmartEdgeBegin(); !it.IsEnd(); ++it) {
			answer[*it] = edge_color_.ColorStr(*it);
		}
		return answer;
	}

	map<VertexId, string> VertexColorMap() const {
		map<VertexId, string> answer;
		for (auto it = this->g().begin(); it != this->g().end(); ++it) {
			answer[*it] = vertex_color_.ColorStr(*it);
		}
		return answer;
	}

	/*virtual*/
	void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
		VERIFY(old_edges.size() > 0);
//		auto color = Color(old_edges.front());
		for (auto it = old_edges.begin(); it != old_edges.end(); ++it) {
//			VERIFY(color == Color(*it));
			Paint(new_edge, Color(*it));
		}
//		Paint(new_edge, color);
	}

	/*virtual*/
	void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		//todo temporary verification
//		VERIFY(Color(edge2) == edge_type::black && new_edge == edge2);
		Paint(new_edge, Color(edge2));
		Paint(new_edge, Color(edge1));
	}

	/*virtual*/
	void HandleSplit(EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge_2) {
		Paint(this->g().EdgeEnd(new_edge_1), Color(old_edge));
		Paint(new_edge_1, Color(old_edge));
		Paint(new_edge_2, Color(old_edge));
	}
};

template<class Graph>
void SaveColoring(const Graph& g
		, const IdTrackHandler<Graph>& int_ids
		, const ColorHandler<Graph>& coloring
		, const string& filename) {
	GraphComponent<Graph> whole_graph(g);
	ofstream stream((filename + ".clr").c_str());
	stream << whole_graph.v_size() << endl;
	for (auto it = whole_graph.v_begin(); it != whole_graph.v_end(); ++it) {
		stream << g.int_id(*it) << " " << int(coloring.Color(*it)) << endl;
	}
	stream << whole_graph.e_size() << endl;
	for (auto it = whole_graph.e_begin(); it != whole_graph.e_end(); ++it) {
		stream << g.int_id(*it) << " " << int(coloring.Color(*it)) << endl;
	}
}

template<class Graph>
void LoadColoring(const Graph& g
		, const IdTrackHandler<Graph>& int_ids
		, ColorHandler<Graph>& coloring
		, const string& filename) {
	ifstream stream((filename + ".clr").c_str());
	size_t v_count;
	stream >> v_count;
	for (size_t i = 0; i < v_count; ++i) {
		size_t id;
		stream >> id;
		size_t color;
		stream >> color;
		coloring.Paint(int_ids.ReturnVertexId(id), color);
	}
	size_t e_count;
	stream >> e_count;
	for (size_t i = 0; i < e_count; ++i) {
		size_t id;
		stream >> id;
		size_t color;
		stream >> color;
		coloring.Paint(int_ids.ReturnEdgeId(id), color);
	}
}

template<class Graph>
auto_ptr<GraphColorer<Graph>> ConstructColorer(
		const ColorHandler<Graph>& coloring) {
	return auto_ptr<GraphColorer<Graph>>(
			new CompositeGraphColorer<Graph>(
					new MapColorer<typename Graph::VertexId>(
							coloring.VertexColorMap()),
					new MapColorer<typename Graph::EdgeId>(
							coloring.EdgeColorMap())));
}

template<class Graph>
auto_ptr<GraphColorer<Graph>> ConstructBorderColorer(const Graph& g,
		const ColorHandler<Graph>& coloring) {
	return auto_ptr<GraphColorer<Graph>>(
			new CompositeGraphColorer<Graph>(
					new BorderVertexColorer<Graph>(g),
					new MapColorer<typename Graph::EdgeId>(
							coloring.EdgeColorMap())));
}

}
