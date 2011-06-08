#ifndef GRAPH_PRINTER_HPP_
#define GRAPH_PRINTER_HPP_

#include <iostream>
#include <fstream>
#include "string"
#include "vector"
#include <map>
#include <sstream>
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

void startSimpleGraphRecord(ostream &out, const string &name);

void endGraphRecord(ostream &out);

void recordParameter(ostream &out, const string &name, const string &value);

string constructCell(const string &label, int border, const string &port);

template<typename tVertex>
void recordVertexId(ostream &out, tVertex id) {
	out << "vertex_" << id;
}

template<typename tVertex>
void recordVertex(ostream &out, Vertex<tVertex> &vertex) {
	recordVertexId(out, vertex.id);
	out << "[";
	recordParameter(out, "label", vertex.label);
	out << ",";
	recordParameter(out, "style", "filled");
	out << ",";
	recordParameter(out, "color", "black");
	out << ",";
	recordParameter(out, "fillcolor", vertex.fillColor);
	out << "]" << endl;
}

template<typename tVertex>
void recordEdge(ostream &out, Edge<tVertex> &edge) {
	recordVertexId(out, edge.from);
	out << "out_";
	out << "->";
	recordVertexId(out, edge.to);
	out << "_in";
	out << "[";
	recordParameter(out, "label", edge.label);
	out << ",";
	recordParameter(out, "color", edge.color);
	out << "]" << endl;
}

template<typename tVertex>
void recordSimpleEdge(ostream &out, Edge<tVertex> &edge) {
	recordVertexId(out, edge.from);
	out << "->";
	recordVertexId(out, edge.to);
	out << "[";
	recordParameter(out, "label", edge.label);
	out << ",";
	recordParameter(out, "color", edge.color);
	out << "]" << endl;
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
string constructNodePairId(tVertex u, tVertex v) {
	stringstream ss;
	if (u == v)
		ss << u;
	else if (u > v)
		ss << v << "_" << u;
	else
		ss << u << "_" << v;
	return ss.str();
}

template<typename tVertex>
string constructNodeId(tVertex v) {
	return constructNodePairId(v, v);
}

template<typename tVertex>
string constructComplexNodeId(string pairId, tVertex v) {
	stringstream ss;
	ss << pairId << ":port_" << v;
	return ss.str();
}

template<typename tVertex>
string vertexIdToString(tVertex v) {
	stringstream ss;
	ss << v;
	return ss.str();
}

template<typename tVertex>
string constructTableEntry(tVertex v, const string &label) {
	stringstream ss;
	ss << "<TR>";
	ss << constructCell("", 0, vertexIdToString(v) + "_in");
	ss << constructCell(label, 0, "");
	ss << constructCell("", 0, vertexIdToString(v) + "out_");
	ss << "</TR>\n";
	return ss.str();
}

template<typename tVertex>
string constructReverceTableEntry(tVertex v, const string &label) {
	stringstream ss;
	ss << "<TR>";
	ss << constructCell("", 0, vertexIdToString(v) + "out_");
	ss << constructCell(label, 0, "");
	ss << constructCell("", 0, vertexIdToString(v) + "_in");
	ss << "</TR>\n";
	return ss.str();
}

template<typename tVertex>
string constructComplexNodeLabel(tVertex v1, const string &label1, tVertex v2,
		const string &label2) {
	return "<TABLE>\n" + constructTableEntry(v1, label1)
			+ constructReverceTableEntry(v2, label2) + "</TABLE>";
}

template<typename tVertex>
string constructComplexNodeLabel(tVertex v, const string &label) {
	return "<TABLE>\n" + constructTableEntry(v, label) + "</TABLE>";
}

template<typename tVertex>
string constructVertexInPairId(tVertex v, tVertex rc) {
	return constructComplexNodeId(constructNodePairId(v, rc), v);
}

string getColor(int currentLength, int approximateLength);

template<typename VertexId>
class GraphPrinter {
protected:
	const string& name_;
	ostream& out_;
public:
	GraphPrinter(const string &name, ostream &out = cout) : name_(name), out_(out) {
	}

	virtual ~GraphPrinter() {
	}

	virtual void open() = 0;

	virtual void close() = 0;

	virtual void AddVertex(VertexId vertexId, const string &label = " ",
			const string &fillColor = "white") = 0;

	virtual void AddEdge(VertexId fromId, VertexId toId, const string &label = " ",
			const string &color = "black") = 0;

};

template<typename VertexId>
class DotGraphPrinter : public GraphPrinter<VertexId> {
private:
	typedef GraphPrinter<VertexId> super;
public:
	DotGraphPrinter(const string &name, ostream &out = cout) : super(name, out) {

	}

	virtual void open() {
		startSimpleGraphRecord(super::out_, super::name_);
	}

	virtual void AddVertex(VertexId vertexId, const string &label,
			const string &fillColor = "white") {
		Vertex<VertexId> v(vertexId, label, fillColor);
		recordVertex<VertexId> (super::out_, v);
	}

	virtual void AddEdge(VertexId fromId, VertexId toId, const string &label,
			const string &color = "black") {
		Edge<VertexId> e(fromId, toId, label, color);
		recordSimpleEdge<VertexId> (super::out_, e);
	}

	virtual void close() {
		endGraphRecord(super::out_);
	}

	virtual ~DotGraphPrinter() {

	}

};

template<typename VertexId>
class PairedGraphPrinter {
private:
	ostream& out_;
	const string& name_;
	map<VertexId, VertexId> vertexMap;
	pair<VertexId, VertexId> currenVertexId;
public:
	PairedGraphPrinter(const string &name, ostream &out = cout)
	: out_(out), name_(name) {
	}

	PairedGraphPrinter(const string &name, const char* filename) {
		out_ = new ofstream(filename, ios::out);
	}

	void open() {
		startGraphRecord(out_, name_);
	}

	void close() {
		endGraphRecord(out_);
	}

	void AddVertex(VertexId v1, string label1, VertexId v2, string label2,
			const string &fillColor = "white") {
		string pairId = constructNodePairId(v1, v2);
		string pairLabel = constructComplexNodeLabel(v1, label1, v2, label2);
		Vertex<string> v(pairId, pairLabel, fillColor);
		recordVertex<string> (out_, v);
	}

	void AddVertex(VertexId v1, string label1, const string &fillColor = "white") {
		string vertexId = constructNodeId(v1);
		string vertexLabel = constructComplexNodeLabel(v1, label1);
		Vertex<string> v(vertexId, vertexLabel, fillColor);
		recordVertex<string> (out_, v);
	}

	void AddEdge(pair<VertexId, VertexId> v1, pair<VertexId, VertexId> v2,
			const string label = " ", const string &color = "black") {
		string v1Id = constructVertexInPairId(v1.first, v1.second);
		string v2Id = constructVertexInPairId(v2.first, v2.second);
		Edge<string> edge(v1Id, v2Id, label, color);
		recordEdge(out_, edge);
	}

	void AddEdge(pair<VertexId, VertexId> v1, VertexId v2,
			const string label = " ", const string &color = "black") {
		string v1Id = constructVertexInPairId(v1.first, v1.second);
		string v2Id = constructVertexInPairId(v2, v2);
		Edge<string> edge(v1Id, v2Id, label, color);
		recordEdge(out_, edge);
	}

	void AddEdge(VertexId v1, pair<VertexId, VertexId> v2,
			const string label = " ", const string &color = "black") {
		string v1Id = constructVertexInPairId(v1, v1);
		string v2Id = constructVertexInPairId(v2.first, v2.second);
		Edge<string> edge(v1Id, v2Id, label, color);
		recordEdge(out_, edge);
	}

	void AddEdge(VertexId v1, VertexId v2, const string label = " ",
			const string &color = "black") {
		string v1Id = constructVertexInPairId(v1, v1);
		string v2Id = constructVertexInPairId(v2, v2);
		Edge<string> edge(v1Id, v2Id, label, color);
		recordEdge(out_, edge);
	}

};


template<class Graph, typename VertexId = typename Graph::VertexId>
class DotPairedGraphPrinter : public GraphPrinter<VertexId> {
private:
	typedef GraphPrinter<VertexId> super;
	const Graph& g_;
	PairedGraphPrinter<VertexId> paired_printer_;
public:
	DotPairedGraphPrinter(const Graph& g, const string &name, ostream &out = cout) : super(name, out), g_(g), paired_printer_(name, out) {

	}

	virtual ~DotPairedGraphPrinter() {

	}

	virtual void open() {
		paired_printer_.open();
	}

	virtual void close() {
		paired_printer_.close();
	}

	virtual void AddVertex(VertexId v, const string &label,
			const string& fillColor) {
		paired_printer_.AddVertex(v, label, g_.conjugate(v), label);
	}

	virtual void AddEdge(VertexId v1, VertexId v2, const string &label,
			const string& color) {
		paired_printer_.AddEdge(make_pair(v1, g_.conjugate(v1)), make_pair(v2, g_.conjugate(v2)), label, color);
	}

};

}
#endif //GRAPH_PRINTER_HPP_//
