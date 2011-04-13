#include "cute.h"
#include <set>
#include "seq.hpp"
#include "sequence.hpp"
#include "nucl.hpp"
#include "read.hpp"
#include "condensed_graph.hpp"
#include "condensed_graph_constructor.hpp"
#include "debruijn.hpp"
#include "graphVisualizer.hpp"
#include <fstream>
#include <tr1/unordered_set>
#include <ext/functional>
#include "test_utils.hpp"

namespace condensed_graph {

using namespace std;

typedef tr1::unordered_set<pair<string, string> , PairHash<string> > edge_set;

typedef tr1::unordered_set<string> vertex_set;

class ToStringHandler: public Traversal::Handler {
	vertex_set& vertices_;
	edge_set& edges_;
public:
	ToStringHandler(vertex_set& vertices, edge_set& edges) :
		vertices_(vertices), edges_(edges) {
	}

	virtual void HandleEdge(const Vertex* v1, const Vertex* v2) {
		edges_.insert(make_pair(v1->nucls().str(), v2->nucls().str()));
	}

	virtual void HandleStartVertex(const Vertex* v) {
		vertices_.insert(v->nucls().str());
	}
};

//void go(const CondensedGraph& g, Vertex* v, set<Vertex*>& visited, string& log) {
//	log += "Entering vertex '" + v->nucls().str() + "'; ";
//	if (visited.count(v) == 1) {
//		log += "Vertex '" + v->nucls().str() + "' has been visited; ";
//	} else {
//		visited.insert(v);
//		vector<Vertex*> desc = g.RightNeighbours(v);
//		for (size_t i = 0; i < desc.size(); ++i) {
//			go(g, desc[i], visited, log);
//		}
//	}
//	log += "Leaving vertex '" + v->nucls().str() + "'; ";
//}

//string printDfs(const CondensedGraph& g, Vertex* start) {
//	DFS dfs(g);
//
//	string log;
//	set<Vertex*> visited;
//	go(g, start, visited, log);
//	return log;
//}

//Vertex* find(const set<Vertex*>& vs, string s) {
//	for (set<Vertex*>::iterator i = vs.begin(); i != vs.end(); ++i) {
//		if ((*i)->nucls().str() == s) {
//			return *i;
//		}
//	}
//	assert(false);
//	return NULL;
//}

//bool contains(const set<Vertex*>& vs, string s) {
//	for (set<Vertex*>::iterator i = vs.begin(); i != vs.end(); ++i) {
//		if ((*i)->nucls().str() == s) {
//			return true;
//		}
//	}
//	return false;
//}

//string print(const set<Vertex*>& vs) {
//	string s = "";
//	for (set<Vertex*>::iterator i = vs.begin(); i != vs.end(); ++i) {
//		s += "'" + (*i)->nucls().str() + "' ";
//	}
//	return s;
//}

string print(const vertex_set& vs) {
	string s = "Vertex set : {";
	for (vertex_set::const_iterator i = vs.begin(); i != vs.end(); ++i) {
		s += "'" + *i + "'; ";
	}
	return s + "}";
}
string print(const edge_set& es) {
	string s = "Edge set : {";
	for (edge_set::const_iterator i = es.begin(); i != es.end(); ++i) {
		s += "'" + (*i).first + "'->'" + (*i).second + "'; ";
	}
	return s;
}

void TestVertex() {
	Vertex* v = new Vertex(Sequence("AAAA"));
	delete v;
}

//void VisTool() {
//	CondensedGraph g(5);
//	Vertex* v1 = g.AddVertex(Sequence("AAAAT"));
//	Vertex* v2 = g.AddVertex(Sequence("AAATA"));
//	Vertex* v3 = g.AddVertex(Sequence("AAATC"));
//	g.LinkVertices(v1, v2);
//	g.LinkVertices(v1, v3);
//	fstream filestr;
//	filestr.open("test.txt", fstream::out);
//	gvis::GraphPrinter<const Vertex*> gp("test graph", filestr);
//	SimpleGraphVisualizer gv(gp);
//	gv.Visualize(g);
//	filestr.close();
//	DFS dfs(g);
//	SimpleStatCounter h;
//	dfs.Traverse(h);
//	cerr << h.v_count() << " " << h.e_count();
//}

//todo refactor
void MyEquals(edge_set es, string s[][2], size_t length) {
	edge_set etalon_edges;
	for (size_t i = 0; i < length; ++i) {
		ASSERT(es.count(make_pair(s[i][0], s[i][1])) == 1);
		ASSERT(es.count(make_pair(complement(s[i][1]), complement(s[i][0]))) == 1);
		etalon_edges.insert(make_pair(s[i][0], s[i][1]));
		etalon_edges.insert(make_pair(complement(s[i][1]), complement(s[i][0])));
	}
	ASSERT_EQUAL(etalon_edges.size(), es.size());
}

//todo refactor
void MyEquals(vertex_set vs, string s[], size_t length) {
	vertex_set etalon_vertices;
	for (size_t i = 0; i < length; ++i) {
		ASSERT(vs.count(s[i]) == 1);
		ASSERT(vs.count(complement(s[i])) == 1);
		etalon_vertices.insert(s[i]);
		etalon_vertices.insert(complement(s[i]));
	}
	ASSERT_EQUAL(etalon_vertices.size(), vs.size());
}

//template <size_t kmer_size_>
//void AssertGraph(size_t read_cnt, string reads_str[], size_t vertex_cnt, string et_vertices[], size_t edge_cnt, string et_edges[][2]) {
//	const vector<Read> reads = MakeReads(reads_str, read_cnt);
//	DirectConstructor<kmer_size_> g_c(reads);
//	CondensedGraph *g;
//	typename DirectConstructor<kmer_size_>::Index *index;
//	g_c.ConstructGraph(g, index);
//
//	edge_set edges;
//	vertex_set vertices;
//	ToStringHandler h(vertices, edges);
//	DFS dfs(*g);
//	dfs.Traverse(h);
//	if (vertex_cnt != 0) {
//		MyEquals(vertices, et_vertices, vertex_cnt);
//	}
//	if (edge_cnt != 0) {
//		MyEquals(edges, et_edges, edge_cnt);
//	}
//
//	delete g;
//	delete index;
//}
//
//template <size_t kmer_size_>
//void AssertGraph(size_t read_cnt, string reads[], size_t vertex_cnt, string et_vertices[]) {
//	AssertGraph<kmer_size_>(read_cnt, reads, vertex_cnt, et_vertices, 0, NULL);
//}
//
//template <size_t kmer_size_>
//void AssertGraph(size_t read_cnt, string reads[], size_t edge_cnt, string et_edges[][2]) {
//	AssertGraph<kmer_size_>(read_cnt, reads, 0, NULL, edge_cnt, et_edges);
//}

template <size_t kmer_size_>
void AssertCondense(size_t read_cnt, string reads_str[], size_t vertex_cnt, string et_vertices[], size_t edge_cnt, string et_edges[][2]) {
	vector<Read> reads = MakeReads(reads_str, read_cnt);
	DeBruijn<kmer_size_> debruijn;
	debruijn.ConstructGraph(reads) ;
	CondenseConstructor<kmer_size_> g_c(debruijn);
	CondensedGraph *g;
	typename CondenseConstructor<kmer_size_>::Index *index;
	g_c.ConstructGraph(g, index);

	edge_set edges;
	vertex_set vertices;
	ToStringHandler h(vertices, edges);
	DFS dfs(*g);
	dfs.Traverse(h);

	DEBUG(print(vertices));
	DEBUG(print(edges));

	if (vertex_cnt != 0) {
		MyEquals(vertices, et_vertices, vertex_cnt);
	}
	if (edge_cnt != 0) {
		MyEquals(edges, et_edges, edge_cnt);
	}

	delete g;
	delete index;
}

template <size_t kmer_size_>
void AssertGraph(size_t read_cnt, string reads[], size_t vertex_cnt, string et_vertices[]) {
	AssertCondense<kmer_size_>(read_cnt, reads, vertex_cnt, et_vertices, 0, NULL);
}

template <size_t kmer_size_>
void AssertGraph(size_t read_cnt, string reads[], size_t edge_cnt, string et_edges[][2]) {
	AssertCondense<kmer_size_>(read_cnt, reads, 0, NULL, edge_cnt, et_edges);
}

void TestSimpleThread() {
	static const size_t read_cnt = 1;
	string reads[read_cnt] = {"ACAAACCACCA"};
	static const size_t vertex_cnt = 1;
	string vertices[vertex_cnt] = {"ACAAACCACCA"};
	AssertGraph<5>(read_cnt, reads, vertex_cnt, vertices);
}

void TestSimpleThread2() {
	static const size_t read_cnt = 2;
	string reads[read_cnt] = {"ACAAACCACCC", "AAACCACCCAC"};
	static const size_t vertex_cnt = 1;
	string vertices[vertex_cnt] = {"ACAAACCACCCAC"};
	AssertGraph<5>(read_cnt, reads, vertex_cnt, vertices);
}

void TestSplitThread() {
	static const size_t read_cnt = 2;
	string reads[read_cnt] = {"ACAAACCACCA", "ACAAACAACCC"};
	static const size_t edge_cnt = 2;
	string edges[edge_cnt][2] = {{"ACAAAC", "AAACAACCC"}, {"ACAAAC", "AAACCACCA"}};
	AssertGraph<5>(read_cnt, reads, edge_cnt, edges);
}

void TestSplitThread2() {
	static const size_t read_cnt = 2;
	string reads[read_cnt] = {"ACAAACCACCA", "ACAAACAACCA"};
	static const size_t edge_cnt = 4;
	string edges[edge_cnt][2] = {{"ACAAAC", "AAACAACC"}, {"AAACAACC", "AACCACCA"}, {"ACAAAC", "AAACC"}, {"AAACC", "AACCACCA"}};
	AssertGraph<5>(read_cnt, reads, edge_cnt, edges);
}

void TestBuldge() {
	static const size_t read_cnt = 2;
	string reads[read_cnt] = {"ACAAAACACCA", "ACAAACCACCA"};
	static const size_t edge_cnt = 4;
	string edges[edge_cnt][2] = {{"ACAAA", "CAAAACACC"}, {"CAAAACACC", "CACCA"}, {"ACAAA", "CAAACCACC"}, {"CAAACCACC", "CACCA"}};
	AssertGraph<5>(read_cnt, reads, edge_cnt, edges);
}

//void TestSimpleHashTable() {
//	SimpleIndex<5> h;
//	Seq<5> k1("AACCG");
//	Vertex* v = new Vertex(Sequence("AAAAAAAAAAAAA"));
//	h.put(k1, v, 1);
//	ASSERT(h.contains(k1));
//	ASSERT_EQUAL(h.get(k1).first, v);
//	ASSERT_EQUAL(h.get(k1).second, 1);
//	h.put(k1, v, 2);
//	ASSERT_EQUAL(h.get(k1).first, v);
//	ASSERT_EQUAL(h.get(k1).second, 2);
//	delete v;
//}

//void TestAddVertex() {
//	//	Graph g;
//	//	g.AddVertex()
//}
//

void TestCondenseSimple() {
	static const size_t read_cnt = 4;

	string reads[read_cnt] = {"CGAAACCAC", "CGAAAACAC", "AACCACACC", "AAACACACC"};
	static const size_t edge_cnt = 4;
	string edges[edge_cnt][2] = {{"CGAAA", "GAAAACACA"}, {"CGAAA", "GAAACCACA"}, {"GAAACCACA", "CACACC"}, {"GAAAACACA", "CACACC"}};

	AssertGraph<5>(read_cnt, reads, edge_cnt, edges);
}

cute::suite CondensedGraphSuite() {
	cute::suite s;
	s.push_back(CUTE(TestVertex));
	 s.push_back(CUTE(TestSimpleThread));
	 s.push_back(CUTE(TestSimpleThread2));
	 s.push_back(CUTE(TestSplitThread));
	 s.push_back(CUTE(TestSplitThread2));
	 s.push_back(CUTE(TestBuldge));
//	 s.push_back(CUTE(TestSimpleHashTable));
	 s.push_back(CUTE(TestCondenseSimple));
//	 s.push_back(CUTE(VisTool));
	return s;
}

}

