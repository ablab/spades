#ifndef GRAPH_VIS_
#define GRAPH_VIS_

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
	out << "_out";
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
	ss << constructCell("", 0, vertexIdToString(v) + "_out");
	ss << "</TR>\n";
	return ss.str();
}

template<typename tVertex>
string constructReverceTableEntry(tVertex v, const string &label) {
	stringstream ss;
	ss << "<TR>";
	ss << constructCell("", 0, vertexIdToString(v) + "_out");
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

template<typename tVertex>
class GraphPrinter {
private:
	ostream *_out;
	map<tVertex, tVertex> vertexMap;
	int approximateLength_;
	int currentLength;
	tVertex currentVertex;
public:
	GraphPrinter(const string &name, ostream &out = cout) {
		_out = &out;
		startSimpleGraphRecord(*_out, name);
	}

	void addVertex(tVertex vertexId, const string &label,
			const string &fillColor = "white") {
		Vertex<tVertex> v(vertexId, label, fillColor);
		recordVertex<tVertex> (*_out, v);
	}

	void addEdge(tVertex fromId, tVertex toId, const string &label = " ",
			const string &color = "black") {
		Edge<tVertex> e(fromId, toId, label, color);
		recordSimpleEdge<tVertex> (*_out, e);
	}

	void output() {
		endGraphRecord(*_out);
	}

	void threadStart(tVertex v, int approximateLength) {
		approximateLength_ = approximateLength;
		currentLength = 0;
		currentVertex = v;
	}

	void threadAdd(tVertex v) {
		addEdge(currentVertex, v, " ",
				getColor(currentLength, approximateLength_));
		currentVertex = v;
		currentLength++;
	}
};

template<typename tVertex>
class PairedGraphPrinter {
private:
	ostream *_out;
	map<tVertex, tVertex> vertexMap;
	int approximateLength_;
	int currentLength;
	pair<tVertex, tVertex> currentVertex;
public:
	PairedGraphPrinter(const string &name, ostream &out = cout) {
		approximateLength_ = -1;
		_out = &out;
		startGraphRecord(*_out, name);
	}

	PairedGraphPrinter(const string &name, const char* filename) {
		_out = new ofstream(filename, ios::out);
		startGraphRecord(*_out, name);
	}

	void addVertex(tVertex v1, string label1, tVertex v2, string label2,
			const string &fillColor = "white") {
		string pairId = constructNodePairId(v1, v2);
		string pairLabel = constructComplexNodeLabel(v1, label1, v2, label2);
		Vertex<string> v(pairId, pairLabel, fillColor);
		recordVertex<string> (*_out, v);
	}

	void addVertex(tVertex v1, string label1, const string &fillColor = "white") {
		string vertexId = constructNodeId(v1);
		string vertexLabel = constructComplexNodeLabel(v1, label1);
		Vertex<string> v(vertexId, vertexLabel, fillColor);
		recordVertex<string> (*_out, v);
	}

	void addEdge(pair<tVertex, tVertex> v1, pair<tVertex, tVertex> v2,
			const string label = " ", const string &color = "black") {
		string v1Id = constructVertexInPairId(v1.first, v1.second);
		string v2Id = constructVertexInPairId(v2.first, v2.second);
		Edge<string> edge(v1Id, v2Id, label, color);
		recordEdge(*_out, edge);
	}

	void addEdge(pair<tVertex, tVertex> v1, tVertex v2,
			const string label = " ", const string &color = "black") {
		string v1Id = constructVertexInPairId(v1.first, v1.second);
		string v2Id = constructVertexInPairId(v2, v2);
		Edge<string> edge(v1Id, v2Id, label, color);
		recordEdge(*_out, edge);
	}

	void addEdge(tVertex v1, pair<tVertex, tVertex> v2,
			const string label = " ", const string &color = "black") {
		string v1Id = constructVertexInPairId(v1, v1);
		string v2Id = constructVertexInPairId(v2.first, v2.second);
		Edge<string> edge(v1Id, v2Id, label, color);
		recordEdge(*_out, edge);
	}

	void addEdge(tVertex v1, tVertex v2, const string label = " ",
			const string &color = "black") {
		string v1Id = constructVertexInPairId(v1, v1);
		string v2Id = constructVertexInPairId(v2, v2);
		Edge<string> edge(v1Id, v2Id, label, color);
		recordEdge(*_out, edge);
	}

	void output() {
		endGraphRecord(*_out);
	}

	void threadStart(pair<tVertex, tVertex> v, int approximateLength) {
		approximateLength_ = approximateLength;
		currentLength = 0;
		currentVertex = v;
	}

	void threadAdd(pair<tVertex, tVertex> v) {
		addEdge(currentVertex, v, " ",
				getColor(currentLength, approximateLength_));
		currentVertex = v;
		currentLength++;
	}
};
}

#endif //GRAPH_VIS_//
