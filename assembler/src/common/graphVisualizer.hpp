#ifndef GRAPH_VIS_
#define GRAPH_VIS_

#include "iostream"
#include "string"
#include "vector"
using namespace std;

namespace gvis {

template<typename tVertex>
struct Vertex {
	tVertex id;
	string label;
	string fillColor;
	Vertex(tVertex _id, string _label, string _fillColor) {
		id = _id;
		label = _label;
		fillColor = _fillColor;
	}
};

template<typename tVertex>
struct Edge {
	tVertex from;
	tVertex to;
	string label;
	string color;
	Edge(tVertex _from, tVertex _to, string _label, string _color) {
		from = _from;
		to = _to;
		label = _label;
		color = _color;
	}
};

void startGraphRecord(ostream &out, const string &name);

void endGraphRecord(ostream &out);

string generateParameterString(const string &name, const string &value);

template<typename tVertex>
void recordVertex(ostream &out, Vertex<tVertex> &vertex) {
	out << vertex.id << " [" << generateParameterString("label", vertex.label)
			<< "," << generateParameterString("fillcolor", vertex.fillColor)
			<< "]" << endl;
}

template<typename tVertex>
void recordEdge(ostream &out, Edge<tVertex> &edge) {
	out << edge.from << "->" << edge.to << "[" << generateParameterString(
			"label", edge.label) << "," << generateParameterString("color",
			edge.color) << "]" << endl;
}

template<typename tVertex>
void recordVertices(ostream &out, vector<Vertex<tVertex> > &vertices) {
	for (typename vector<Vertex<tVertex> >::iterator it = vertices.begin(); it
			!= vertices.end(); it++) {
		recordVertex(out, *it);
	}
}

template<typename tVertex>
void recordEdges(ostream &out, vector<Edge<tVertex> > &edges) {
	for (typename vector<Edge<tVertex> >::iterator it = edges.begin(); it
			!= edges.end(); it++) {
		recordEdge(out, *it);
	}
}

template<typename tVertex>
void outputGraph(ostream &out, const string &graphName,
		vector<Vertex<tVertex> > &vertices, vector<Edge<tVertex> > &edges) {
	startGraphRecord(out, graphName);
	recordVertices<tVertex> (out, vertices);
	recordEdges<tVertex> (out, edges);
	endGraphRecord(out);
}

template<typename tVertex>
class IGraphPrinter {
protected:
	ostream *_out;
public:
	virtual void addVertex(tVertex vertexId, const string &label,
			const string &fillColor = "white") = 0;
	virtual void addEdge(tVertex fromId, tVertex toId, const string &label,
			const string &color = "black") = 0;
	virtual void output() = 0;
};

template<typename tVertex>
class GraphPrinter: public IGraphPrinter<tVertex> {
public:
	GraphPrinter(const string &name, ostream &out) {
		IGraphPrinter<tVertex>::_out = &out;
		startGraphRecord(*IGraphPrinter<tVertex>::_out, name);
	}

	GraphPrinter(const string &name) {
		IGraphPrinter<tVertex>::_out = &cout;
		startGraphRecord(*IGraphPrinter<tVertex>::_out, name);
	}

	virtual void addVertex(tVertex vertexId, const string &label,
			const string &fillColor = "white") {
		Vertex<tVertex> v(vertexId, label, fillColor);
		recordVertex<tVertex> (*IGraphPrinter<tVertex>::_out, v);
	}

	virtual void addEdge(tVertex fromId, tVertex toId, const string &label,
			const string &color = "black") {
		Edge<tVertex> e(fromId, toId, label, color);
		recordEdge<tVertex> (*IGraphPrinter<tVertex>::_out, e);
	}

	void output() {
		endGraphRecord(*IGraphPrinter<tVertex>::_out);
	}
};

//template<typename tVertex>
//class GraphScheme: public IGraphPrinter<tVertex> {
//private:
//	string _name;
//	vector<Vertex<tVertex> > _vertices;
//	vector<Edge<tVertex> > _edges;
//public:
//
//	GraphScheme(string name) {
//		_name = name;
//		IGraphPrinter<tVertex>::_out = &cout;
//	}
//
//	GraphScheme(string name, ostream &out) {
//		_name = name;
//		IGraphPrinter<tVertex>::_out = &out;
//	}
//
//	virtual void addVertex(tVertex vertexId, const string &label,
//			const string &fillColor = "white") {
//		Vertex<tVertex> v(vertexId, label, fillColor);
//		_vertices.push_back(v);
//	}
//
//	virtual void addEdge(tVertex fromId, tVertex toId, const string &label,
//			const string &color = "black") {
//		Edge<tVertex> e(fromId, toId, label, color);
//		_edges.push_back(e);
//	}
//
//	virtual void output() {
//		outputGraph<tVertex> (*IGraphPrinter<tVertex>::_out, _name, _vertices,
//				_edges);
//	}
//
//	void output(ostream out) {
//		outputGraph<tVertex> (out, _name, _vertices, _edges);
//	}
//};
}

#endif //GRAPH_VIS_//
