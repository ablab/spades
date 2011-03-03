#ifndef GRAPH_VIS_
#define GRAPH_VIS_

#include "iostream"
#include "string"
#include "vector"
using namespace std;

namespace gvis {

void startGraphRecord(ostream &out, string &name);

void endGraphRecord(ostream &out);

template<typename tVertex>
void recordVertex(ostream &out, tVertex vertexId, string &vertexLabel) {
	if (vertexLabel.length() != 0)
		out << vertexId << " [label" << "=" << vertexLabel << "]" << endl;
	else
		out << vertexId << ";" << endl;
}

template<typename tVertex>
void recordEdge(ostream &out, tVertex fromId, tVertex toId, string &edgeLabel) {
	out << fromId << "->" << toId;
	if (edgeLabel.length() != 0)
		out << "[label=" << edgeLabel << "]";
	out << endl;
}

template<typename tVertex>
void recordVertices(ostream &out, vector<pair<tVertex, string> > &vertices) {
	for (typename vector<pair<tVertex, string> >::iterator it =
			vertices.begin(); it != vertices.end(); it++) {
		pair<tVertex, string> v = *it;
		recordVertex(out, v.first, v.second);
	}
}

template<typename tVertex>
void recordEdges(ostream &out,
		vector<pair<pair<tVertex, tVertex> , string> > &edges) {
	for (typename vector<pair<pair<tVertex, tVertex> , string> >::iterator it =
			edges.begin(); it != edges.end(); it++) {
		pair<pair<tVertex, tVertex> , string> e = *it;
		recordEdge(out, e.first.first, e.first.second, e.second);
	}
}

template<typename tVertex>
void outputGraph(ostream &out, string &graphName,
		vector<pair<tVertex, string> > &vertices,
		vector<pair<pair<tVertex, tVertex> , string> > &edges) {
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
	void addVertex(tVertex vertexId, const string &label);
	void addEdge(tVertex fromId, tVertex toId, const string &label);
	void output();
};

template<typename tVertex>
class OnlineGraphPrinter: public IGraphPrinter<tVertex> {
	OnlineGraphPrinter(string name, ostream &out) {
		IGraphPrinter<tVertex>::_out = &out;
		startGraphRecord(IGraphPrinter<tVertex>::_out, name);
	}

	OnlineGraphPrinter(string name) {
		IGraphPrinter<tVertex>::_out = &cout;
		startGraphRecord(IGraphPrinter<tVertex>::_out, name);
	}

	void addVertex(tVertex vertexId, const string &label) {
		recordVertex(IGraphPrinter<tVertex>::_out, vertexId, label);
	}

	void addEdge(tVertex fromId, tVertex toId, const string &label) {
		recordEdge(IGraphPrinter<tVertex>::_out, fromId, toId, label);
	}

	void output() {
		endGraphRecord(IGraphPrinter<tVertex>::_out);
	}
};

template<typename tVertex>
class GraphScheme: public IGraphPrinter<tVertex> {
private:
	string _name;
	vector<pair<tVertex, string> > _vertices;
	vector<pair<pair<tVertex, tVertex> , string> > _edges;
public:

	GraphScheme(string name) {
		_name = name;
		IGraphPrinter<tVertex>::_out = &cout;
	}

	GraphScheme(string name, ostream &out) {
		_name = name;
		IGraphPrinter<tVertex>::_out = &out;
	}

	void addVertex(tVertex vertexId, const string &label) {
		_vertices.push_back(make_pair(vertexId, label));
	}

	void addEdge(tVertex fromId, tVertex toId, const string &label) {
		_edges.push_back(make_pair(make_pair(fromId, toId), label));
	}

	void output() {
		outputGraph<tVertex> (*IGraphPrinter<tVertex>::_out, _name, _vertices,
				_edges);
	}

	void output(ostream out) {
		outputGraph<tVertex> (out, _name, _vertices, _edges);
	}
};
}

#endif //GRAPH_VIS_//
