#ifndef EDGEGRAPHTEST_HPP_
#define EDGEGRAPHTEST_HPP_
#include "edge_graph_constructor.hpp"
#include "test_utils.hpp"
#include "cute.h"
#include "strobe_reader.hpp"
#include "paired_info.hpp"
#include "simple_tools.hpp"
#include "debruijn_plus.hpp"
#include <tr1/unordered_set>
//LOGGER("d.edge_graph_test");

namespace edge_graph {

using de_bruijn::Traversal;
using de_bruijn::DFS;
using de_bruijn::PairedInfoIndex;

void EmptyGraphTest() {
	EdgeGraph g(11);
	ASSERT_EQUAL(11, g.k());
	ASSERT_EQUAL(0u, g.size());
}

void OneVertexGraphTest() {
	EdgeGraph g(11);
	g.AddVertex();
	ASSERT_EQUAL(2u, g.size());
	Vertex *v = *(g.begin());
	Vertex *rcv = g.Complement(v);
	ASSERT(v != rcv);
	ASSERT_EQUAL(v, g.Complement(rcv));
}

pair<vector<VertexId> , vector<EdgeId> > createGraph(EdgeGraph &graph,
		int edgeNumber) {
	vector<VertexId> v;
	vector<EdgeId> e;
	v.push_back(graph.AddVertex());
	for (int i = 0; i < edgeNumber; i++) {
		v.push_back(graph.AddVertex());
		e.push_back(
				graph.AddEdge(v[v.size() - 2], v[v.size() - 1],
						Sequence("AAAAAAAAAAAAAAAAA")));
	}
	return make_pair(v, e);
}

void OneEdgeGraphTest() {
	EdgeGraph g(11);
	pair<vector<VertexId> , vector<EdgeId> > data = createGraph(g, 1);
	ASSERT_EQUAL(1u, g.OutgoingEdgeCount(data.first[0]));
	ASSERT_EQUAL(0u, g.OutgoingEdgeCount(data.first[1]));
	ASSERT_EQUAL(data.second[0], g.GetUniqueOutgoingEdge(data.first[0]));
	ASSERT_EQUAL(g.Complement(data.second[0]),
			g.GetUniqueOutgoingEdge(g.Complement(data.first[1])));
	ASSERT_EQUAL(data.second[0],
			g.Complement(g.Complement(data.second[0])));
	ASSERT_EQUAL(!(g.EdgeNucls(data.second[0])),
			g.EdgeNucls(g.Complement(data.second[0])));
}

void EdgeMethodsSimpleTest() {
	EdgeGraph g(11);
	pair<vector<VertexId> , vector<EdgeId> > data = createGraph(g, 2);
	ASSERT_EQUAL(data.second[0], &g.GetData(data.second[0]));
	ASSERT_EQUAL(
			true,
			g.AreLinkable(data.first[0], data.first[1],
					Sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")));
	ASSERT_EQUAL(
			false,
			g.AreLinkable(data.first[0], data.first[1],
					Sequence("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")));
}

void VertexMethodsSimpleTest() {
	EdgeGraph g(11);
	pair<vector<VertexId> , vector<EdgeId> > data = createGraph(g, 2);
	ASSERT_EQUAL(data.second[0], g.GetUniqueIncomingEdge(data.first[1]));
	ASSERT_EQUAL(data.second[0], g.GetUniqueOutgoingEdge(data.first[0]));
	ASSERT_EQUAL(false, g.CanCompressVertex(data.first[0]));
	ASSERT_EQUAL(true, g.CanCompressVertex(data.first[1]));
	ASSERT_EQUAL(false, g.CheckUniqueIncomingEdge(data.first[0]));
	ASSERT_EQUAL(true, g.CheckUniqueIncomingEdge(data.first[1]));
	ASSERT_EQUAL(false, g.CheckUniqueOutgiongEdge(data.first[2]));
	ASSERT_EQUAL(true, g.CheckUniqueOutgiongEdge(data.first[1]));
	ASSERT_EQUAL(true, g.IsDeadEnd(data.first[2]));
	ASSERT_EQUAL(false, g.IsDeadEnd(data.first[1]));
	ASSERT_EQUAL(true, g.IsDeadStart(data.first[0]));
	ASSERT_EQUAL(false, g.IsDeadStart(data.first[1]));
}

//void GraphMethodsSimpleTest() {
//	EdgeGraph g(11);
//	pair<vector<VertexId> , vector<EdgeId> > data = createGraph(g, 2);
//	ASSERT_EQUAL(vector<ActionHandler*> (), g.GetHandlers());
//	ActionHandler* handler = new ActionHandler();
//	g.AddActionHandler(handler);
//	vector<ActionHandler*> handlers = g.GetHandlers();
//	ASSERT_EQUAL(1u, handlers.size());
//	ASSERT_EQUAL(handler, handlers[0]);
//	g.RemoveActionHandler(handler);
//	ASSERT_EQUAL(vector<ActionHandler*> (), g.GetHandlers());
//}

void SmartIteratorTest() {
	EdgeGraph g(11);
	pair<vector<VertexId> , vector<EdgeId> > data = createGraph(g, 4);
	size_t num = 0;
	set<VertexId> visited;
	std::less<VertexId> comp;
	SmartVertexIterator<EdgeGraph> it = g.SmartVertexBegin(comp);
	SmartVertexIterator<EdgeGraph> it1 = g.SmartVertexBegin(comp);
	SmartVertexIterator<EdgeGraph> it2 = g.SmartVertexEnd(comp);
	SmartVertexIterator<EdgeGraph> it3 = g.SmartVertexEnd(comp);
	for (SmartVertexIterator<EdgeGraph> it = g.SmartVertexBegin(comp); g.SmartVertexEnd(
			comp) != it; ++it) {
		num++;
		DEBUG( "with seq in vert" << g.VertexNucls(*it).str());
		visited.insert(*it);
	}
	ASSERT_EQUAL(num, data.first.size() * 2);
	for (size_t i = 0; i < data.first.size(); i++) {
		ASSERT(visited.find(data.first[i]) != visited.end());
		ASSERT(visited.find(g.Complement(data.first[i])) != visited.end());
	}
}

typedef de_bruijn::PairedInfoIndex<EdgeGraph>::PairInfo PairInfo;
typedef string MyRead;
typedef pair<MyRead, MyRead> MyPairedRead;
typedef string MyEdge;
typedef pair<MyEdge, MyEdge> MyEdgePair;
typedef map<MyEdgePair, pair<size_t, double>> EdgePairInfo;
typedef map<MyEdge, double> CoverageInfo;
typedef tr1::unordered_set<MyEdge> Edges;

string print(const Edges& es) {
	string s = "Edge set : {";
	for (auto i = es.begin(); i != es.end(); ++i) {
		s += "'" + *i + "'; ";
	}
	return s;
}

class ToStringHandler: public TraversalHandler {
	Edges& edges_;
	EdgeGraph g_;
public:
	ToStringHandler(Edges& edges, EdgeGraph &g) :
		edges_(edges), g_(g) {
	}

	virtual void HandleEdge(EdgeId e) {
		//todo rewrite using graph object (maybe add g_ to superclass)
		edges_.insert(g_.EdgeNucls(e).str());
	}

};

const Edges AddComplement(const Edges& edges) {
	Edges ans;
	for (auto it = edges.begin(); it != edges.end(); ++it) {
		ans.insert(*it);
		ans.insert(Complement(*it));
	}
	return ans;
}

const CoverageInfo AddComplement(const CoverageInfo& coverage_info) {
	CoverageInfo ans;
	for (auto it = coverage_info.begin(); it != coverage_info.end(); ++it) {
		ans.insert(*it);
		ans.insert(make_pair(Complement((*it).first), (*it).second));
	}
	return ans;
}

void EdgesEqual(const Edges& s1, const Edges& s2) {
	ASSERT_EQUAL(s1.size(), s2.size());
	for (auto it = s1.begin(); it != s1.end(); ++it) {
		ASSERT(s2.count(*it) > 0);
	}
}

const vector<Read> MakeReads(const vector<MyRead>& reads) {
	vector<Read> ans;
	for (size_t i = 0; i < reads.size(); ++i) {
		ans.push_back(Read("", reads[i], ""));
	}
	return ans;
}

const vector<PairedRead> MakePairedReads(const vector<MyPairedRead>& paired_reads, size_t insert_size) {
	vector<PairedRead> ans;
	for (size_t i = 0; i < paired_reads.size(); ++i) {
		ans.push_back(PairedRead(Read("", paired_reads[i].first, ""), Read("", paired_reads[i].second, ""), insert_size));
	}
	return ans;
}

void AssertEdges(EdgeGraph& g, const Edges& etalon_edges) {
	Edges edges;
	for (auto it = g.SmartEdgeBegin(); it != g.SmartEdgeEnd(); ++it) {
		edges.insert(g.EdgeNucls(*it).str());
	}
	EdgesEqual(edges, etalon_edges);
}

template<size_t kmer_size_>
void AssertGraph(const vector<string>& reads, const vector<string>& etalon_edges) {
	typedef VectorStream<Read> Stream;
	Stream read_stream(MakeReads(reads));
	EdgeGraph g(kmer_size_);
	de_bruijn::EdgeIndex<kmer_size_ + 1, EdgeGraph> index(g);

	ConstructGraph<kmer_size_, Stream>(g, index, read_stream);

	AssertEdges(g, AddComplement(Edges(etalon_edges.begin(), etalon_edges.end())));
}

bool EqualDouble(double d1, double d2) {
	return std::abs(d1 - d2) < 1e-5;
}

void AssertCoverage(EdgeGraph& g, const CoverageInfo& etalon_coverage) {
	for (auto it = g.SmartEdgeBegin(); it != g.SmartEdgeEnd(); ++it) {
		CoverageInfo::const_iterator cov_info_it = etalon_coverage.find(g.EdgeNucls(*it).str());
		ASSERT(cov_info_it != etalon_coverage.end() && EqualDouble(g.coverage(*it), (*cov_info_it).second));
	}
}

void AssertPairInfo(const PairedIndex& paired_index, const EdgePairInfo& etalon_pair_info) {

}

template<size_t k>
void AssertGraph(const vector<MyPairedRead>& paired_reads, size_t insert_size, const vector<MyEdge>& etalon_edges
		, const CoverageInfo& etalon_coverage, const EdgePairInfo& etalon_pair_info) {
	typedef VectorStream<PairedRead> PairedReadStream;

	PairedReadStream paired_read_stream(MakePairedReads(paired_reads, insert_size));
	EdgeGraph g(k);
	EdgeIndex<k + 1, EdgeGraph> index(g);
	de_bruijn::CoverageHandler<EdgeGraph> coverage_handler(g);
	PairedIndex paired_index(g);

	ConstructGraphWithPairedInfo<k, PairedReadStream>(g, index, coverage_handler, paired_index, paired_read_stream);

	AssertEdges(g, AddComplement(Edges(etalon_edges.begin(), etalon_edges.end())));

	AssertCoverage(g, etalon_coverage);

	AssertPairInfo(paired_index, etalon_pair_info);
}

//todo rename tests

void TestSimpleThread() {
	vector<string> reads = { "ACAAACCACCA" };
//	vector<string> edges = { "ACAAACCACCA" };
	AssertGraph<5> (reads, reads);
}

void TestSimpleThread2() {
	vector<string> reads = { "ACAAACCACCC", "AAACCACCCAC" };
	vector<string> edges = { "ACAAACCACCCAC" };
	AssertGraph<5> (reads, edges);
}

void TestSplitThread() {
	vector<string> reads = { "ACAAACCACCA", "ACAAACAACCC" };
	vector<string> edges = { "ACAAAC", "CAAACCACCA", "CAAACAACCC" };
	AssertGraph<5> (reads, edges);
}

void TestSplitThread2() {
	vector<string> reads = { "ACAAACCACCA", "ACAAACAACCA" };
	vector<string> edges = { "AACCACCA", "ACAAAC", "CAAACCA", "CAAACAACCA" };
	AssertGraph<5> (reads, edges);
}

void TestBuldge() {
	vector<string> reads = { "ACAAAACACCA", "ACAAACCACCA" };
//	vector<string> edges = { "ACAAAACACCA", "ACAAACCACCA" };
	AssertGraph<5> (reads, reads);
}

void TestCondenseSimple() {
	vector<string> reads = { "CGAAACCAC", "CGAAAACAC", "AACCACACC", "AAACACACC" };
	vector<string> edges = { "CGAAAACACAC", "CACACC", "CGAAACCACAC" };
	AssertGraph<5> (reads, edges);
}

void TestPairedInfo() {
	vector<string> reads = { "CGAAACCAC", "CGAAAACAC", "AACCACACC", "AAACACACC" };
	vector<string> edges = { "CGAAAACACAC", "CACACC", "CGAAACCACAC" };
	AssertGraph<5> (reads, edges);
}

cute::suite EdgeGraphSuite() {
	cute::suite s;
	s.push_back(CUTE(EmptyGraphTest));
	s.push_back(CUTE(OneVertexGraphTest));
	s.push_back(CUTE(OneEdgeGraphTest));
	s.push_back(CUTE(EdgeMethodsSimpleTest));
	s.push_back(CUTE(VertexMethodsSimpleTest));
//	s.push_back(CUTE(GraphMethodsSimpleTest));
	s.push_back(CUTE(SmartIteratorTest));
	s.push_back(CUTE(TestBuldge));

	s.push_back(CUTE(TestSimpleThread));
	s.push_back(CUTE(TestSimpleThread2));
	s.push_back(CUTE(TestSplitThread));
	s.push_back(CUTE(TestSplitThread2));
	s.push_back(CUTE(TestCondenseSimple));
	return s;
}
}

#endif /* EDGEGRAPHTEST_HPP_ */
