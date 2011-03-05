#include "cute.h"
#include <set>
#include "seq.hpp"
#include "sequence.hpp"
#include "nucl.hpp"
#include "condensedGraph.hpp"

namespace condensed_graph {

using namespace std;

bool contains(const set<Vertex*>& vs, string s) {
	for (set<Vertex*>::iterator i = vs.begin(); i != vs.end(); ++i) {
		if ((*i)->nucls().str() == s) {
			return true;
		}
	}
	return false;
}

string print(const set<Vertex*>& vs) {
	string s = "";
	for (set<Vertex*>::iterator i = vs.begin(); i != vs.end(); ++i) {
		s += "'" + (*i)->nucls().str() + "' " ;
	}
	return s;
}

void TestVertex() {
	Vertex* v = new Vertex(Sequence("AAAA"));
	delete v;
}

void TestSimpleThread() {
	Graph g;
	Read r("ACAAACCACC");//"ATGCATATGC");
	DEBUG("Read is " + r.str())
	g.ThreadRead(r);
	//repeat for complication
	g.ThreadRead(r);
	ASSERT(g.component_roots().size() == 2);
	ASSERT(contains(g.component_roots(), r.str()));
	ASSERT(contains(g.component_roots(), (!r).str()));
}

std::string complement(const std::string& s) {
	return (!Sequence(s)).str();
}

void TestSplitThread() {
	Graph g;
	Read r1("ACAAACCACC");//"ATGCATATGC");
	Read r2("ACAAACAACC");//"ATGCATATGC");
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
	ASSERT(contains(g.component_roots(), complement("AAACCACC")));
	ASSERT(contains(g.component_roots(), complement("AAACAACC")));
}

void TestSimpleHashTable() {
	SimpleHashTable h;
	Kmer k1("AACCG");
	Vertex* null_v = (Vertex*) NULL;
	h.put(k1, make_pair(null_v, 1));
	ASSERT(h.contains(k1));
	ASSERT_EQUAL(h.get(k1).first, null_v);
	ASSERT_EQUAL(h.get(k1).second, 1);

	h.put(k1, make_pair(null_v, 2));
	ASSERT_EQUAL(h.get(k1).first, null_v);
	ASSERT_EQUAL(h.get(k1).second, 2);
}

}

using namespace condensed_graph;
cute::suite CondensedGraphSuite() {
	cute::suite s;
	s.push_back(CUTE(TestVertex));
	s.push_back(CUTE(TestSimpleHashTable));
	s.push_back(CUTE(TestSimpleThread));
	s.push_back(CUTE(TestSplitThread));
	return s;
}
