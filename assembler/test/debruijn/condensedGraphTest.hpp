#include "cute.h"
#include <set>
#include "seq.hpp"
#include "sequence.hpp"
#include "nucl.hpp"
#include "condensedGraph.hpp"

namespace condensed_graph {

using namespace std;

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

void TestSimpleThread() {
	Graph g;
	Read r("ACAAACCACCA");//"ATGCATATGC");
	DEBUG("Read is " + r.str())
	g.ThreadRead(r);
	//repeat for complication
	g.ThreadRead(r);
	DEBUG(print(g.component_roots()))
	ASSERT(g.component_roots().size() == 2);
	ASSERT(contains(g.component_roots(), r.str()));
	ASSERT(contains(g.component_roots(), (!r).str()));
	ASSERT_EQUAL(printDfs(g, find(g.component_roots(), "ACAAACCACCA"))
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

	DEBUG(print(g.component_roots()))
	ASSERT(g.component_roots().size() == 2);
	ASSERT(contains(g.component_roots(), "ACAAACCACCCAC"));
	ASSERT(contains(g.component_roots(), complement("ACAAACCACCCAC")));
	ASSERT_EQUAL(printDfs(g, find(g.component_roots(), "ACAAACCACCCAC"))
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
	DEBUG(print(g.component_roots()))
	ASSERT(g.component_roots().size() == 3);
	ASSERT(contains(g.component_roots(), "ACAAAC"));
	ASSERT(contains(g.component_roots(), complement("AAACCACCA")));
	ASSERT(contains(g.component_roots(), complement("AAACAACCC")));
	ASSERT_EQUAL(printDfs(g, find(g.component_roots(), "ACAAAC"))
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
	DEBUG(print(g.component_roots()))
	ASSERT(g.component_roots().size() == 2);
	ASSERT(contains(g.component_roots(), "ACAAAC"));
	ASSERT(contains(g.component_roots(), complement("AACCACCA")));
	ASSERT_EQUAL(printDfs(g, find(g.component_roots(), "ACAAAC"))
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
	//repeat for complication
	g.ThreadRead(r1);
	g.ThreadRead(r2);
	DEBUG(print(g.component_roots()))
	ASSERT(g.component_roots().size() == 2);
	ASSERT(contains(g.component_roots(), "ACAAA"));
	ASSERT(contains(g.component_roots(), complement("CACCA")));
	ASSERT_EQUAL(printDfs(g, find(g.component_roots(), "ACAAA"))
			, "Entering vertex 'ACAAA'; Entering vertex 'CAAAACACC'; Entering vertex 'CACCA'; Leaving vertex 'CACCA'; Leaving vertex 'CAAAACACC'; Entering vertex 'CAAACCACC'; Entering vertex 'CACCA'; Vertex 'CACCA' has been visited; Leaving vertex 'CACCA'; Leaving vertex 'CAAACCACC'; Leaving vertex 'ACAAA'; "
	);
}

void TestSimpleHashTable() {
	SimpleHashTable h;
	Kmer k1("AACCGAACCGA");
	Vertex* v = new Vertex(Sequence("AAAA"));
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

}

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
	return s;
}
