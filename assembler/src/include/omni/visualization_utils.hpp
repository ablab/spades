#ifndef VISUALIZATIONUTILS_HPP_
#define VISUALIZATIONUTILS_HPP_

#include "graph_printer.hpp"
#include "graph_labeler.hpp"
#include "omni_utils.hpp"
#include <stack>
#include <queue>
#include "dijkstra.hpp"
#include "splitters.hpp"
#include "abstract_conjugate_graph.hpp"
#include "abstract_nonconjugate_graph.hpp"

namespace omnigraph {

//DECL_LOGGER("omg.gvis")

using omnigraph::GraphLabeler;
using omnigraph::EmptyGraphLabeler;
using omnigraph::SmartEdgeIterator;
using omnigraph::Path;
using omnigraph::UnorientedDijkstra;

template<class Graph>
class GraphVisualizer {
	typedef typename Graph::VertexId VertexId;
protected:
	const Graph& g_;
	//	const GraphLabeler<Graph> &labeler_;
	//	gvis::GraphPrinter<VertexId>& gp_;
public:
	GraphVisualizer(const Graph& g/*, gvis::GraphPrinter<VertexId>& gp*/) :
			g_(g) {

	}

	virtual ~GraphVisualizer() {

	}

	virtual void Visualize() = 0;

};

template<class Graph>
class PartialGraphVisualizer {
	typedef typename Graph::VertexId VertexId;
protected:
	const Graph& g_;
	GraphPrinter<VertexId>& gp_;
public:
	PartialGraphVisualizer(const Graph& g, GraphPrinter<VertexId>& gp) :
			g_(g), gp_(gp) {
	}

	void open() {
		gp_.open();
	}

	void close() {
		gp_.close();
	}

	virtual ~PartialGraphVisualizer() {

	}

	virtual void Visualize(const vector<VertexId>& vertices) = 0;

};

template<class Graph>
class SimpleGraphVisualizer: public GraphVisualizer<Graph> {
	typedef GraphVisualizer<Graph> super;
	typedef typename Graph::VertexId VertexId;
	GraphPrinter<VertexId>& gp_;
	const omnigraph::GraphLabeler<Graph>& gl_;
public:
	SimpleGraphVisualizer(const Graph& g, GraphPrinter<VertexId>& gp,
			const omnigraph::GraphLabeler<Graph>& gl) :
			super(g), gp_(gp), gl_(gl) {
	}

	virtual void Visualize() {
		gp_.open();
		TRACE("Visualize started");
		for (auto it = super::g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
			gp_.AddVertex(*it, gl_.label(*it));
		}TRACE("Vertices printed");
		for (auto it = super::g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			gp_.AddEdge(super::g_.EdgeStart(*it), super::g_.EdgeEnd(*it),
					gl_.label(*it), "black", super::g_.length(*it));
		}TRACE("Edges printed");
		gp_.close();
	}

private:
	DECL_LOGGER("SimpleGraphVisualizer")
};

template<class Graph>
class AdapterGraphVisualizer: public GraphVisualizer<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef GraphVisualizer<Graph> super;
	PartialGraphVisualizer<Graph>& partial_visualizer_;
public:
	AdapterGraphVisualizer(const Graph& g,
			PartialGraphVisualizer<Graph>& partial_visualizer) :
			super(g), partial_visualizer_(partial_visualizer) {

	}

	virtual void Visualize() {
		partial_visualizer_.open();
		partial_visualizer_.Visualize(
				vector<VertexId>(super::g_.begin(), super::g_.end()));
		partial_visualizer_.close();
	}
};

template<class Graph, typename ElementId>
class ElementColorer {
public:
	virtual string GetColour(ElementId element) const = 0;

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

	const auto_ptr<ElementColorer<Graph, VertexId>> vertex_colorer_;
	const auto_ptr<ElementColorer<Graph, EdgeId>> edge_colorer_;
public:
//	CompositeGraphColorer(auto_ptr<ElementColorer<Graph, VertexId>> vertex_colorer
//			, auto_ptr<ElementColorer<Graph, EdgeId>> edge_colorer) :
//				vertex_colorer_(vertex_colorer),
//				edge_colorer_(edge_colorer) {
//
//	}

	CompositeGraphColorer(ElementColorer<Graph, VertexId>* vertex_colorer
			, ElementColorer<Graph, EdgeId>* edge_colorer) :
				vertex_colorer_(vertex_colorer),
				edge_colorer_(edge_colorer) {

	}

	/*virtual */string GetColour(VertexId v) const {
		return vertex_colorer_->GetColour(v);
	}

	/*virtual */map<VertexId, string> GetColours(const set<VertexId> &vertices) const {
		return vertex_colorer_->GetColours(vertices);
	}

	/*virtual */string GetColour(EdgeId e) const {
		return edge_colorer_->GetColour(e);
	}

	/*virtual */map<EdgeId, string> GetColours(const set<EdgeId> &edges) const {
		return edge_colorer_->GetColours(edges);
	}

};

template<class Graph, typename ElementId>
class MapColorer : public ElementColorer<Graph, ElementId> {
private:
	map<ElementId, string> color_map_;
	optional<string> default_color_;
public:
	MapColorer(const map<ElementId, string> &color_map) : color_map_(color_map) {
	}

	virtual ~MapColorer() {
	}

	MapColorer(const map<ElementId, string> &color_map, const string& default_color) :
		color_map_(color_map),
		default_color_(default_color) {
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
class BorderVertexColorer : public ElementColorer<Graph, typename Graph::VertexId> {
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
		VERIFY(false);
		return "";
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
class ColoredGraphVisualizer: public PartialGraphVisualizer<Graph> {
	typedef PartialGraphVisualizer<Graph> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const omnigraph::GraphLabeler<Graph>& gl_;
	const GraphColorer<Graph> &colorer_;

public:
	ColoredGraphVisualizer(const Graph& g, GraphPrinter<VertexId>& gp,
			const omnigraph::GraphLabeler<Graph>& gl,
			const GraphColorer<Graph> &colorer) :
			super(g, gp), gl_(gl), colorer_(colorer) {
	}

	virtual void Visualize(const vector<VertexId>& vertices) {
		set<VertexId> vertex_set(vertices.begin(), vertices.end());
		map<VertexId, string> vertex_colour_map = colorer_.GetColours(vertex_set);

		for (auto v_it = vertex_set.begin(); v_it != vertex_set.end(); ++v_it) {
			super::gp_.AddVertex(*v_it, gl_.label(*v_it), vertex_colour_map[*v_it]);
		}
		for (auto v_it = vertex_set.begin(); v_it != vertex_set.end(); ++v_it) {
			const vector<EdgeId> edges = super::g_.OutgoingEdges(*v_it);
			TRACE("Working with vertex " << *v_it);
			for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it) {
				VertexId edge_end = super::g_.EdgeEnd(*e_it);
				TRACE(
						super::g_.coverage(*e_it) << " " << super::g_.length(*e_it));
				if (vertex_set.count(edge_end) > 0) {
					super::gp_.AddEdge(*v_it, edge_end, gl_.label(*e_it),
							colorer_.GetColour(*e_it), super::g_.length(*e_it));
					TRACE("Edge added");
				}
			}
		}
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

// edge_colorer management is passed here
template <class Graph>
auto_ptr<GraphColorer<Graph>> DefaultColorer(const Graph& g,
		ElementColorer<Graph, typename Graph::EdgeId>* edge_colorer) {
	return auto_ptr<GraphColorer<Graph>>(new CompositeGraphColorer<Graph>(new BorderVertexColorer<Graph>(g), edge_colorer));
}

template <class Graph>
auto_ptr<GraphColorer<Graph>> DefaultColorer(const Graph& g,
		const map<typename Graph::EdgeId, string>& edge_color_map,
		const string& default_color = "") {
	return DefaultColorer(g, new MapColorer<Graph, typename Graph::EdgeId>(edge_color_map, default_color));
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
void WriteToDotFile(const Graph& g, const GraphLabeler<Graph>& labeler,
		const string& file_name, const string& graph_name = "my_graph") {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	DotGraphPrinter<typename Graph::VertexId> gpr(graph_name, filestr);
	SimpleGraphVisualizer<Graph> sgv(g, gpr, labeler);
	sgv.Visualize();
	filestr.close();
}

template<class Graph>
void WriteSimple(const Graph& g, const GraphLabeler<Graph>& labeler,
		const string& file_name, const string& graph_name) {
	fstream filestr;
	string simple_file_name(file_name);
	//	simple_file_name.insert(simple_file_name.size() - 4, "_simple");
	filestr.open((simple_file_name).c_str(), fstream::out);
	DotGraphPrinter<typename Graph::VertexId> gpr(graph_name, filestr);
	SimpleGraphVisualizer<Graph> sgv(g, gpr, labeler);
	sgv.Visualize();
	filestr.close();
}

template<class Graph>
void WriteSimple(const Graph& g, const GraphLabeler<Graph>& labeler,
		const string& file_name, const string& graph_name,
		const Path<typename Graph::EdgeId> &path1,
		const Path<typename Graph::EdgeId> &path2) {
	fstream filestr;
	string simple_file_name(file_name);
	//	simple_file_name.insert(simple_file_name.size() - 4, "_simple");
	filestr.open(simple_file_name.c_str(), fstream::out);
	DotGraphPrinter<typename Graph::VertexId> gp(graph_name, filestr);
	auto_ptr<GraphColorer<Graph>> colorer(DefaultColorer(g, path1, path2));
	ColoredGraphVisualizer<Graph> gv(g, gp, labeler, *colorer);
	AdapterGraphVisualizer<Graph> result_vis(g, gv);
	result_vis.Visualize();
	filestr.close();
}

template<class Graph>
void WritePaired(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& file_name,
		const string& graph_name,
		const Path<typename Graph::EdgeId> &path1/* = Path<typename Graph::EdgeId> ()*/,
		const Path<typename Graph::EdgeId> &path2/* = Path<typename Graph::EdgeId> ()*/) {
	typedef typename Graph::EdgeId EdgeId;
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	DotPairedGraphPrinter<Graph> gp(g, graph_name, labeler, filestr);
	CompositeGraphColorer<Graph> colorer(/*create_auto_ptr(*/new BorderVertexColorer<Graph>(g)
			, /*create_auto_ptr(*/new MapColorer<Graph, EdgeId>(PathColorer<Graph>(g, path1, path2).ColorPath(), ""));
	ColoredGraphVisualizer<Graph> gv(g, gp, labeler, colorer);
	AdapterGraphVisualizer<Graph> result_vis(g, gv);
	result_vis.Visualize();
	filestr.close();
}

template<class Graph>
class VisualizerFactory {
public:
	virtual auto_ptr<PartialGraphVisualizer<Graph>> GetVisualizerInstance(
			GraphPrinter<typename Graph::VertexId> &gp) = 0;
	virtual auto_ptr<GraphPrinter<typename Graph::VertexId>> GetPrinterInstance(
			const string &graph_name, ostream &os) = 0;
	virtual ~VisualizerFactory() {
	}
};

template<class Graph>
class ColoredVisualizerFactory: public VisualizerFactory<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	const Graph &graph_;
	const GraphLabeler<Graph> &labeler_;
	const GraphColorer<Graph> &colorer_;

	auto_ptr<GraphPrinter<VertexId>> PrinterInstance(
			const Graph& graph, const string &graph_name,
			ostream &os) {
		return auto_ptr<GraphPrinter<VertexId>>(new DotPairedGraphPrinter<Graph>(
				graph, graph_name, labeler_, os));
	}

	auto_ptr<GraphPrinter<VertexId>> PrinterInstance(
			const AbstractNonconjugateGraph<typename Graph::DataMaster>& graph, const string &graph_name,
			ostream &os, const GraphLabeler<Graph> &labeler) {
		return auto_ptr<GraphPrinter<VertexId>>(new DotGraphPrinter<VertexId>(graph_name, os));
	}

public:
	ColoredVisualizerFactory(const Graph& graph,
			const GraphLabeler<Graph> &labeler,
			const GraphColorer<Graph> &colorer) :
			graph_(graph), labeler_(labeler), colorer_(colorer) {
	}

	virtual auto_ptr<GraphPrinter<typename Graph::VertexId>> GetPrinterInstance(
			const string &graph_name, ostream &os) {
		return PrinterInstance(graph_, graph_name, os);
	}

	virtual auto_ptr<PartialGraphVisualizer<Graph>> GetVisualizerInstance(
			GraphPrinter<typename Graph::VertexId> &gp) {
		return auto_ptr<PartialGraphVisualizer<Graph>>(new ColoredGraphVisualizer<Graph>(graph_, gp, labeler_,
				colorer_));
	}

	virtual ~ColoredVisualizerFactory() {
	}
};

template<class Graph>
class ComponentGraphVisualizer: public GraphVisualizer<Graph> {
private:
	VisualizerFactory<Graph> &factory_;
	ComponentSplitter<typename Graph::VertexId> &splitter_;
	const string &file_name_;
	const string &graph_name_;
	size_t max_parts_number_;

	string ConstructComponentName(const string &file_name, size_t cnt,
			const string &component_name) {
		stringstream ss;
		ss << "_" << cnt << component_name;
		string res = file_name;
		res.insert(res.length() - 4, ss.str());
		return res;
	}

public:
	ComponentGraphVisualizer(const Graph &graph,
			VisualizerFactory<Graph> &factory,
			ComponentSplitter<typename Graph::VertexId> &splitter,
			const string &file_name, const string &graph_name,
			size_t max_parts_number = 100) :
			GraphVisualizer<Graph>(graph), factory_(factory), splitter_(
					splitter), file_name_(file_name), graph_name_(graph_name), max_parts_number_(
					max_parts_number) {
	}

	virtual ~ComponentGraphVisualizer() {
	}

	virtual void Visualize() {
		size_t cnt = 1;
		while (!splitter_.Finished() && cnt <= max_parts_number_) {
			auto component = splitter_.NextComponent();
			string component_name = ConstructComponentName(file_name_, cnt,
					splitter_.ComponentName());
			ofstream os;
			os.open(component_name.c_str());
			auto_ptr<GraphPrinter<typename Graph::VertexId>> gp =
					factory_.GetPrinterInstance(graph_name_, os);
			auto_ptr<PartialGraphVisualizer<Graph>> visualizer = factory_.GetVisualizerInstance(*gp);
			visualizer->open();
			if (component.size() >= 1000) { // what the magic constant???
				WARN("Too large component " << component.size());
			}
			visualizer->Visualize(component);
			visualizer->close();
			os.close();
			cnt++;
		}
	}
};

template<class Graph>
string InsertComponentName(string file_name, string component) {
	string res = file_name;
	res.insert(res.length() - 4, "_" + component + "_");
	return res;
}

template<class Graph>
void WriteErrors(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& file_name,
		const string& graph_name,
		const Path<typename Graph::EdgeId> &path1/* = Path<typename Graph::EdgeId> ()*/,
		const Path<typename Graph::EdgeId> &path2/* = Path<typename Graph::EdgeId> ()*/) {
	PathColorer<Graph> path_colorer(g, path1, path2);
	set<typename Graph::EdgeId> black = path_colorer.BlackEdges();
	ErrorComponentSplitter<Graph> splitter(g, black);
	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
	ComponentSizeFilter<Graph> checker(g, 1500, 2);
	string error_file_name = InsertComponentName<Graph>(file_name, "error");
	WriteComponents(g, splitter, checker, error_file_name, *DefaultColorer(g, coloring), labeler);

//	PathColorer<Graph> path_colorer(g, path1, path2);
//	set<typename Graph::EdgeId> black = path_colorer.BlackEdges();
//	ErrorComponentSplitter<Graph> splitter(g, black);
//	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
//	ColoredVisualizerFactory<Graph> factory(g, labeler, coloring);
//	string error_file_name = InsertComponentName<Graph>(file_name, "error");
//	ComponentGraphVisualizer<Graph> gv(g, factory, splitter, error_file_name,
//			graph_name, 200);
//	gv.Visualize();
	//	size_t cnt = 0;
	//	while (!splitter.Finished() && cnt < 100) {
	//		fstream filestr;
	//		filestr.open(ConstructComponentName<Graph> (file_name, cnt).c_str(),
	//				fstream::out);
	//		gvis::DotPairedGraphPrinter<Graph> gp(g, graph_name, filestr);
	//		ColoredGraphVisualizer<Graph> gv(g, gp, gl, coloring);
	//		auto component = splitter.NextComponent();
	//		gp.open();
	//		gv.Visualize(component);
	//		gp.close();
	//		cnt++;
	//	}
}

//todo strange similar looking methods!!!
template<class Graph>
void WriteComponents(const Graph& g,
		ComponentSplitter<typename Graph::VertexId> &splitter,
		const string& file_name,
		const GraphColorer<Graph> &colorer,
		const GraphLabeler<Graph>& labeler,
		const string& graph_name = "my_graph") {
	ColoredVisualizerFactory<Graph> factory(g, labeler, colorer);
	ComponentGraphVisualizer<Graph> gv(g, factory, splitter, file_name,
			graph_name, 24000);
	gv.Visualize();
}

template<class Graph>
void WriteComponents(const Graph& g,
		ComponentSplitter<typename Graph::VertexId> &inner_splitter,
		const AbstractFilter<vector<typename Graph::VertexId>> &checker,
		const string& file_name,
		const GraphColorer<Graph> &colorer,
		const GraphLabeler<Graph>& labeler,
		const string& graph_name = "my_graph") {
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, checker);
	WriteComponents<Graph>(g, splitter, file_name, colorer,
			labeler, graph_name);
}

template<class Graph>
void WriteComponents(const Graph& g, size_t split_edge_length,
		const string& file_name,
		const GraphColorer<Graph>& colorer,
		const GraphLabeler<Graph>& labeler,
		const string& graph_name = "my_graph") {
	ReliableSplitter<Graph> splitter(g, 60, split_edge_length);
	ComponentSizeFilter<Graph> filter(g, split_edge_length, 2);
	WriteComponents<Graph>(g, splitter, filter, file_name,
			colorer, labeler, graph_name);
}

//todo alert!!! magic constants!!!
//todo refactoring of params needed
template<class Graph>
void WriteComponentsAlongPath(const Graph& g,
		const GraphLabeler<Graph>& labeler, const string& file_name,
		size_t split_edge_length,
		const MappingPath<typename Graph::EdgeId>& path,
		Path<typename Graph::EdgeId> color1 = Path<typename Graph::EdgeId>(),
		Path<typename Graph::EdgeId> color2 = Path<typename Graph::EdgeId>(), bool colour_path = false) {
//	Path<typename Graph::EdgeId> simple_path1 = color1.simple_path();
//	Path<typename Graph::EdgeId> simple_path2 = color2.simple_path();
	PathColorer<Graph> path_colorer(g, /*simple_path1*/color1, /*simple_path2*/
			color2);
	auto coloring = path_colorer.ColorPath();
	if(colour_path) {
		for(size_t i = 0; i < path.size(); i++) {
			coloring.insert(make_pair(path[i].first, "green"));
		}
	}
	//	LongEdgesSplitter<Graph> inner_splitter(g, split_edge_length);
	//	ReliableSplitterAlongGenome(g, 60, split_edge_length, MappingPath<EdgeId> genome_path)
	ReliableSplitterAlongPath<Graph> splitter(g, 60, split_edge_length,
			path);
	ComponentSizeFilter<Graph> filter(g, 1000000, 0);
	WriteComponents<Graph>(g, splitter, filter, file_name,
			*DefaultColorer(g, coloring), labeler);
}

//todo alert!!! magic constants!!!
template<class Graph>
void WriteComponentsAlongGenome(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& file_name,
		size_t split_edge_length,
		MappingPath<typename Graph::EdgeId> color1 = MappingPath<
				typename Graph::EdgeId>(),
		MappingPath<typename Graph::EdgeId> color2 = MappingPath<
				typename Graph::EdgeId>()) {
	WriteComponentsAlongPath<Graph>(g, labeler, file_name,
			split_edge_length, color1, color1.simple_path(),
			color2.simple_path());
}

template<class Graph>
void WriteComponentsAroundEdge(const Graph& g
		, typename Graph::EdgeId e
		, const string& file_name
		, const GraphColorer<Graph>& colorer
		, const GraphLabeler<Graph>& labeler) {
	EdgeNeighborhoodFinder<Graph> splitter(g, e, 50, 500);
	WriteComponents(g, splitter/*, "locality_of_edge_" + ToString(g_.int_id(edge))*/
			, file_name
			, colorer, labeler);
}

template<class Graph>
void WriteToFile(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& file_name,
		const string& graph_name,
		const Path<typename Graph::EdgeId> &path1/* = Path<typename Graph::EdgeId> ()*/,
		const Path<typename Graph::EdgeId> &path2/* = Path<typename Graph::EdgeId> ()*/) {
//	if (g.size() < 10000) {
	WritePaired(g, labeler, file_name, graph_name, path1, path2);
	WriteSimple(g, labeler, file_name, graph_name, path1, path2);
//	}
	WriteErrors(g, labeler, file_name, graph_name, path1, path2);
}

}

#endif /* VISUALIZATIONUTILS_HPP_ */
