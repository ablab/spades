#include "cute.h"
#include <set>
#include "seq.hpp"
#include "sequence.hpp"
#include "nucl.hpp"
#include "condensedGraph.hpp"
#include "debruijn.hpp"
#include "graphVisualizer.hpp"
#include <fstream>
#include <tr1/unordered_set>
#include <ext/functional>

namespace condensed_graph {

using namespace std;

template <typename T>
struct PairHash {
	size_t operator() (pair<T, T> p) const {
		return hash<T>()(p.first) + hash<T>()(p.second);
	}
};

template <typename T>
struct PairLess {
	bool operator() (pair<T, T> p1, pair<T, T> p2) const {
		return less<T>()(p1.first, p2.first) ? true : (less<T>()(p2.first, p1.first) ? false : less<T>()(p1.second, p2.second));
	}
};

typedef tr1::unordered_set<pair<string, string>, PairHash<string> > edge_set;

class EdgeStringHandler: public Traversal::Handler {
	edge_set& set_;
public:
	EdgeStringHandler(edge_set& set) : set_(set) {

	}
	virtual void HandleEdge(const Vertex* v1, const Vertex* v2) {
		set_.insert(make_pair(v1->nucls().str(), v2->nucls().str()));
	}
};

void go(const Graph& g, Vertex* v, set<Vertex*>& visited, string& log) {
	log += "Entering vertex '" + v->nucls().str() + "'; ";
	if (visited.count(v) == 1) {
		log += "Vertex '" + v->nucls().str() + "' has been visited; ";
	} else {
		visited.insert(v);
		vector<Vertex*> desc = g.Desc(v);
		for (size_t i = 0; i < desc.size(); ++i) {
			go(g, desc[i], visited, log);
		}
	}
	log += "Leaving vertex '" + v->nucls().str() + "'; ";
}

string printDfs(const Graph g, Vertex* start) {
	DFS dfs(g);

	string log;
	set<Vertex*> visited;
	go(g, start, visited, log);
	return log;
}

Vertex* find(const set<Vertex*>& vs, string s) {
	for (set<Vertex*>::iterator i = vs.begin(); i != vs.end(); ++i) {
		if ((*i)->nucls().str() == s) {
			return *i;
		}
	}
	assert(false);
	return NULL;
}

bool contains(const set<Vertex*>& vs, string s) {
	for (set<Vertex*>::iterator i = vs.begin(); i != vs.end(); ++i) {
		if ((*i)->nucls().str() == s) {
			return true;
		}
	}
	return false;
}

std::string complement(const std::string& s) {
	return (!Sequence(s)).str();
}

string print(const set<Vertex*>& vs) {
	string s = "";
	for (set<Vertex*>::iterator i = vs.begin(); i != vs.end(); ++i) {
		s += "'" + (*i)->nucls().str() + "' ";
	}
	return s;
}

void TestVertex() {
	Vertex* v = new Vertex(Sequence("AAAA"));
	delete v;
}

void VisTool() {
	Graph g;
	g.AddVertex(Sequence("AAAAT"));
	g.AddVertex(Sequence("AAATA"));
	g.AddVertex(Sequence("AAATC"));
	g.LinkVertices(g.GetPosition(Kmer("AAAAT")).first,
			g.GetPosition(Kmer("AAATA")).first);
	g.LinkVertices(g.GetPosition(Kmer("AAAAT")).first,
			g.GetPosition(Kmer("AAATC")).first);
	fstream filestr;
	filestr.open("test.txt", fstream::out);
	gvis::GraphPrinter<const Vertex*> gp("test graph", filestr);
	SimpleGraphVisualizer gv(gp);
	gv.Visualize(g);
	filestr.close();
	DFS dfs(g);
	SimpleStatCounter h;
	dfs.Traverse(h);
	cerr<<h.v_count()<<" "<<h.e_count();
}

void TestSimpleThread() {
	Graph g;
	Read r("ACAAACCACCA");//"ATGCATATGC");
	DEBUG("Read is " + r.str())
	g.ThreadRead(r);
	//repeat for complication
	g.ThreadRead(r);
	DEBUG(print(g.vertices()))
	ASSERT(g.vertices().size() == 2);
	ASSERT(contains(g.vertices(), r.str()));
	ASSERT(contains(g.vertices(), (!r).str()));
	ASSERT_EQUAL(printDfs(g, find(g.vertices(), "ACAAACCACCA"))
			,"Entering vertex 'ACAAACCACCA'; Leaving vertex 'ACAAACCACCA'; "
	);
}

void TestSimpleThread2() {
	Graph g;
	Read r1("ACAAACCACCC");
	Read r2("AAACCACCCAC");
	DEBUG("Read 1 is " + r1.str())
	g.ThreadRead(r1);
	//repeat for complication
	g.ThreadRead(r1);

	DEBUG("Read 2 is " + r2.str())
	g.ThreadRead(r2);
	//repeat for complication
	g.ThreadRead(r2);

	DEBUG(print(g.vertices()))
	ASSERT(g.vertices().size() == 2);
	ASSERT(contains(g.vertices(), "ACAAACCACCCAC"));
	ASSERT(contains(g.vertices(), complement("ACAAACCACCCAC")));
	ASSERT_EQUAL(printDfs(g, find(g.vertices(), "ACAAACCACCCAC"))
			,"Entering vertex 'ACAAACCACCCAC'; Leaving vertex 'ACAAACCACCCAC'; "
	);
}

void TestSplitThread() {
	//for n = 11 AACCA - repeat
	Graph g;
	Read r1("ACAAACCACCA");
	Read r2("ACAAACAACCC");
	DEBUG("Read 1 is " + r1.str())
	DEBUG("Read 2 is " + r2.str())
	g.ThreadRead(r1);
	g.ThreadRead(r2);
	//repeat for complication
	g.ThreadRead(r1);
	g.ThreadRead(r2);
	DEBUG(print(g.vertices()));
	ASSERT(g.vertices().size() == 6);
	ASSERT(contains(g.vertices(), "ACAAAC"));
	ASSERT(contains(g.vertices(), complement("AAACCACCA")));
	ASSERT(contains(g.vertices(), complement("AAACAACCC")));
	ASSERT_EQUAL(printDfs(g, find(g.vertices(), "ACAAAC"))
			,"Entering vertex 'ACAAAC'; Entering vertex 'AAACAACCC'; Leaving vertex 'AAACAACCC'; Entering vertex 'AAACCACCA'; Leaving vertex 'AAACCACCA'; Leaving vertex 'ACAAAC'; "
	);
}

void TestSplitThread2() {
	//for n = 11 AACCA - repeat
	Graph g;
	Read r1("ACAAACCACCA");
	Read r2("ACAAACAACCA");
	DEBUG("Read 1 is " + r1.str())
	DEBUG("Read 2 is " + r2.str())
	g.ThreadRead(r1);
	g.ThreadRead(r2);
	//repeat for complication
	g.ThreadRead(r1);
	g.ThreadRead(r2);
	DEBUG(print(g.vertices()));
	ASSERT(g.vertices().size() == 8);
	ASSERT(contains(g.vertices(), "ACAAAC"));
	ASSERT(contains(g.vertices(), complement("AACCACCA")));
	ASSERT_EQUAL(printDfs(g, find(g.vertices(), "ACAAAC"))
			, "Entering vertex 'ACAAAC'; Entering vertex 'AAACAACC'; Entering vertex 'AACCACCA'; Leaving vertex 'AACCACCA'; Leaving vertex 'AAACAACC'; Entering vertex 'AAACC'; Entering vertex 'AACCACCA'; Vertex 'AACCACCA' has been visited; Leaving vertex 'AACCACCA'; Leaving vertex 'AAACC'; Leaving vertex 'ACAAAC'; "
	);
}

void TestBuldge() {
	//for n = 11 AACCA - repeat
	Graph g;
	Read r1("ACAAAACACCA");
	Read r2("ACAAACCACCA");
	DEBUG("Read 1 is " + r1.str())
	DEBUG("Read 2 is " + r2.str())
	g.ThreadRead(r1);
	g.ThreadRead(r2);
	//repeat for complicationTraverse
	g.ThreadRead(r1);
	g.ThreadRead(r2);
	DEBUG(print(g.vertices()))
	ASSERT(g.vertices().size() == 8);
	ASSERT(contains(g.vertices(), "ACAAA"));
	ASSERT(contains(g.vertices(), complement("CACCA")));
	ASSERT_EQUAL(printDfs(g, find(g.vertices(), "ACAAA"))
			, "Entering vertex 'ACAAA'; Entering vertex 'CAAAACACC'; Entering vertex 'CACCA'; Leaving vertex 'CACCA'; Leaving vertex 'CAAAACACC'; Entering vertex 'CAAACCACC'; Entering vertex 'CACCA'; Vertex 'CACCA' has been visited; Leaving vertex 'CACCA'; Leaving vertex 'CAAACCACC'; Leaving vertex 'ACAAA'; "
	);
}

void TestSimpleHashTable() {
	SimpleHashTable h;
	Kmer k1("AACCG");
	Vertex* v = new Vertex(Sequence("AAAAAAAAAAAAA"));
	h.put(k1, v, 1);
	ASSERT(h.contains(k1));
	ASSERT_EQUAL(h.get(k1).first, v);
	ASSERT_EQUAL(h.get(k1).second, 1);
	h.put(k1, v, 2);
	ASSERT_EQUAL(h.get(k1).first, v);
	ASSERT_EQUAL(h.get(k1).second, 2);
	delete v;
}

void TestAddVertex() {
	//	Graph g;
	//	g.AddVertex()
}

void MyEquals(edge_set e, string s[][2], size_t length) {
	set<pair<string, string>, PairLess<string> > etalon_edges;
	for (size_t i = 0; i < length; ++i) {
		ASSERT(e.count(make_pair(s[i][0], s[i][1])) == 1);
		ASSERT(e.count(make_pair(complement(s[i][1]), complement(s[i][0]))) == 1);
		etalon_edges.insert(make_pair(s[i][0], s[i][1]));
		etalon_edges.insert(make_pair(complement(s[i][1]), complement(s[i][0])));
	}
	cout << etalon_edges.size() << endl;
	ASSERT_EQUAL(etalon_edges.size(), e.size());
}

void TestCondenseSimple() {
	string ss[] = {"CGAAACCAC", "CGAAAACAC", "AACCACACC", "AAACACACC"};
	vector<strobe_read<R, 4> > input;
	input.push_back(strobe_read<R, 4>(ss));
	DeBruijn<K> g;
	g.ConstructGraph(input);
	condensed_graph::Graph condensed;
	CondenseGraph(g, condensed);
	edge_set set;
	EdgeStringHandler h(set);
	DFS dfs(condensed);
	dfs.Traverse(h);
	string s[][2] = {{"CGAAA", "GAAAACACA"}, {"CGAAA", "GAAAACACA"}, {"GAAACCACA", "CACACC"}, {"GAAAACACA", "CACACC"}};
	MyEquals(set, s, 4);

	for (edge_set::iterator it = set.begin(); it != set.end(); it++) {
		cout << (*it).first << "  " << (*it).second << endl;
	}
	cout << set.size() << endl;
}

}
/*
using namespace condensed_graph;
cute::suite CondensedGraphSuite() {
	cute::suite s;
	s.push_back(CUTE(TestVertex));
	s.push_back(CUTE(TestSimpleHashTable));
	s.push_back(CUTE(TestSimpleThread));
	s.push_back(CUTE(TestSimpleThread2));
	s.push_back(CUTE(TestBuldge));
	s.push_back(CUTE(TestSplitThread));
	s.push_back(CUTE(TestSplitThread2));
//	s.push_back(CUTE(VisTool));
	return s;
}*/

using namespace condensed_graph;
cute::suite CondensedGraphSuite() {
	cute::suite s;
	s.push_back(CUTE(TestCondenseSimple));
	return s;
}

