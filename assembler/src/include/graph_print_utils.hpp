//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef GRAPH_PRINTER_HPP_
#define GRAPH_PRINTER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>

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
	int len;
	Edge(tVertex _from, tVertex _to, string _label, string _color, int _len) {
		from = _from;
		to = _to;
		label = _label;
		color = _color;
		len = _len;
	}
};

void startGraphRecord(ostream &out, const string &name);

void startSimpleGraphRecord(ostream &out, const string &name);

void endGraphRecord(ostream &out);

void recordParameter(ostream &out, const string &name, const string &value);
void recordParameterInQuotes(ostream &out, const string &name, const string &value); 

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
	recordParameterInQuotes(out, "label", edge.label);
	out << ",";
	recordParameter(out, "len", ToString(edge.len));
	out << ",";
	recordParameter(out, "K", ToString(edge.len));
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
	recordParameterInQuotes(out, "label", edge.label);
	out << ",";
	recordParameter(out, "len", ToString(edge.len));
	out << ",";
	recordParameter(out, "K", ToString(edge.len));
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
string IdToStr(tVertex u) {
	stringstream ss;
	ss << u;
	return ss.str();
}

template<typename tVertex>
string constructNodePairId(tVertex u, tVertex v) {
	stringstream ss;
	string u_str = IdToStr(u);
	string v_str = IdToStr(v);
	if (u == v)
		ss << u;
	else if (u_str > v_str)
		ss << v_str << "_" << u_str;
	else
		ss << u_str << "_" << v_str;
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
string constructTableEntry(tVertex v, const string &label) {
	stringstream ss;
	ss << "<TR>";
	ss << constructCell("", 0, ToString(v) + "_in");
	ss << constructCell(label, 0, "");
	ss << constructCell("", 0, ToString(v) + "_out");
	ss << "</TR>\n";
	return ss.str();
}

template<typename tVertex>
string constructReverceTableEntry(tVertex v, const string &label) {
	stringstream ss;
	ss << "<TR>";
	ss << constructCell("", 0, ToString(v) + "_out");
	ss << constructCell(label, 0, "");
	ss << constructCell("", 0, ToString(v) + "_in");
	ss << "</TR>\n";
	return ss.str();
}

template<typename tVertex>
string constructComplexNodeLabel(tVertex v1, const string &label1, tVertex v2,
		const string &label2, const string &fillColor) {
	return "<TABLE bgcolor = \"" + fillColor + "\">\n" + constructTableEntry(v1, label1)
			+ constructReverceTableEntry(v2, label2) + "</TABLE>";
}

template<typename tVertex>
string constructComplexNodeLabel(tVertex v, const string &label, const string &fillColor) {
	return "<TABLE bgcolor = \"" + fillColor + "\">\n" + constructTableEntry(v, label) + "</TABLE>";
}

template<typename tVertex>
string constructVertexInPairId(tVertex v, tVertex rc) {
	return constructComplexNodeId(constructNodePairId(v, rc), v);
}

string getColor(int currentLength, int approximateLength);

template<typename VertexId>
class PairedGraphPrinter {
private:
	ostream& out_;
	const string name_;
	map<VertexId, VertexId> vertexMap;
	pair<VertexId, VertexId> currenVertexId;
public:
	PairedGraphPrinter(const string &name, ostream &out = cout) :
		out_(out), name_(name) {
	}

	void open() {
		startGraphRecord(out_, name_);
	}

	void close() {
		endGraphRecord(out_);
	}

	void AddVertex(VertexId v1, string label1, VertexId v2, string label2,
			const string fillColor = "white") {
		string pairId = constructNodePairId(v1, v2);
		string pairLabel = constructComplexNodeLabel(v1, label1, v2, label2, fillColor);
		Vertex<string> v(pairId, pairLabel, fillColor);
		recordVertex<string> (out_, v);
	}

	void AddVertex(VertexId v1, string label1,
			const string fillColor = "white") {
		string vertexId = constructNodeId(v1);
		string vertexLabel = constructComplexNodeLabel(v1, label1);
		Vertex<string> v(vertexId, vertexLabel, fillColor);
		recordVertex<string> (out_, v);
	}

	void AddEdge(pair<VertexId, VertexId> v1, pair<VertexId, VertexId> v2,
			const string label = " ", const string color = "black", const int length = 0) {
		string v1Id = constructVertexInPairId(v1.first, v1.second);
		string v2Id = constructVertexInPairId(v2.first, v2.second);
		Edge<string> edge(v1Id, v2Id, label, color, length);
		recordEdge(out_, edge);
	}

	void AddEdge(pair<VertexId, VertexId> v1, VertexId v2,
			const string label = " ", const string color = "black", const int length = 0) {
		string v1Id = constructVertexInPairId(v1.first, v1.second);
		string v2Id = constructVertexInPairId(v2, v2);
		Edge<string> edge(v1Id, v2Id, label, color, length);
		recordEdge(out_, edge);
	}

	void AddEdge(VertexId v1, pair<VertexId, VertexId> v2,
			const string label = " ", const string color = "black", const int length = 0) {
		string v1Id = constructVertexInPairId(v1, v1);
		string v2Id = constructVertexInPairId(v2.first, v2.second);
		Edge<string> edge(v1Id, v2Id, label, color, length);
		recordEdge(out_, edge);
	}

	void AddEdge(VertexId v1, VertexId v2, const string label = " ",
			const string color = "black", const int length = 0) {
		string v1Id = constructVertexInPairId(v1, v1);
		string v2Id = constructVertexInPairId(v2, v2);
		Edge<string> edge(v1Id, v2Id, label, color, length);
		recordEdge(out_, edge);
	}

};

// moved from cpp
inline void startGraphRecord(ostream &out, const string &name) {
	out << "digraph " << name << " {" << endl;
	out << "node" << "[";
	recordParameter(out, "fontname", "Courier");
	out << ",";
	recordParameter(out, "shape", "plaintext");
	out << "]" << endl;
}

inline void startSimpleGraphRecord(ostream &out, const string &name) {
	out << "digraph " << name << " {" << endl;
	out << "node" << "[";
	recordParameter(out, "fontname", "Courier");
	out << "]" << endl;
/* 
 * If this does not break any principles of
 * drawing pictures in debrujin, please uncomment
 * and commit. Multicolored pics are very bad w/o
 * penwidth change
 * fondarat@gmail.com
 *
  out << "edge" << "[";
  recordParameter(out, "penwidth", "1.8");
  out << "]" << endl;
*/
}

inline void endGraphRecord(ostream &out) {
	out << "}" << endl;
}

inline void recordParameter(ostream &out, const string &name, const string &value) {
	out << name << "=" << "<" << value << ">";
}

inline void recordParameterInQuotes(ostream &out, const string &name, const string &value) {
	out << name << "=" << "\"" << value << "\"";
}

inline string constructCell(const string &label, int border, const string &port) {
	stringstream ss;
	ss << "<TD BORDER = \"" << border << "\" PORT = \"port_" << port << "\">"
			<< label << "</TD>";
	return ss.str();
}

inline double getColorParameter(int l, int r, double perc) {
	return l * perc + r * (1 - perc);
}

inline string getColor(int currentLength, int approximateLength) {
	currentLength %= approximateLength;
	int points[8][3] = {{0, 0, 1}, {0, 1, 1}, {1, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 1}, {0, 0, 1}};
	stringstream ss;
	int bound = approximateLength / 6;
	int num = currentLength / bound;
	double perc = (currentLength % bound) * 1. / bound;
	for(int i = 0; i < 3; i++) {
		ss << getColorParameter(points[num][i], points[num + 1][i], perc);
		if(i != 2)
			ss << ",";
	}
	return ss.str();
}

}
#endif //GRAPH_PRINTER_HPP_//
