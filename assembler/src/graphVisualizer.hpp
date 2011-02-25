#ifndef GRAPH_VIS_
#define GRAPH_VIS_

#include "iostream"
#include "string"
#include "vector"
using namespace std;

namespace gvis {

void startGraphRecord(string name);

void endGraphRecord();

template<typename tVertex>
void recordVertex(tVertex vertexId, string vertexLabel) {
	if(vertexLabel.length() != 0)
		cout << vertexId << " [label" << "=" << vertexLabel << "]" << endl;
	else
		cout << vertexId << ";" << endl;
}

template<typename tVertex>
void recordEdge(tVertex fromId, tVertex toId, string edgeLabel) {
	cout << fromId << "->" << toId;
	if(edgeLabel.length() != 0)
		cout << "[label=" << edgeLabel << "]";
	cout << endl;
}

template<typename tVertex>
void recordVertices(vector<pair<tVertex, string> > vertices) {
	for(typename vector<pair<tVertex, string> >::iterator it = vertices.begin(); it != vertices.end(); it++) {
		pair<tVertex, string> v = *it;
		recordVertex(v.first, v.second);
	}
}

template<typename tVertex>
void recordEdges(vector<pair<pair<tVertex, tVertex>, string> > edges) {
	for(typename vector<pair<pair<tVertex, tVertex>, string> >::iterator it = edges.begin(); it != edges.end(); it++) {
		pair<pair<tVertex, tVertex>, string> e = *it;
		recordEdge(e.first.first, e.first.second, e.second);
	}
}

template<typename tVertex>
void outputGraph(string graphName, vector<pair<tVertex, string> > vertices, vector< pair< pair<tVertex, tVertex>, string> > edges) {
	startGraphRecord(graphName);
	recordVertices<tVertex>(vertices);
	recordEdges<tVertex>(edges);
	endGraphRecord();
}

template<typename tVertex>
class GraphScheme {
public:
	string _name;
	vector<pair<tVertex, string> > _vertices;
	vector< pair< pair<tVertex, tVertex>, string> > _edges;

	GraphScheme(string name) {
		_name = name;
	}

	void addVertex(tVertex vertexId, string label) {
		_vertices.push_back(make_pair(vertexId, label));
	}

	void addEdge(tVertex fromId, tVertex toId, string label) {
		_edges.push_back(make_pair(make_pair(fromId, toId), label));
	}

	void output() {
		outputGraph<tVertex>(_name, _vertices, _edges);
	}
};
}
#endif //GRAPH_VIS_//
