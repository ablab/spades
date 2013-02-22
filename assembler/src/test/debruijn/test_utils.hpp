//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <boost/test/unit_test.hpp>
#include "launch.hpp"
#include "graph_construction.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/vector_reader.hpp"
#include "io/converting_reader_wrapper.hpp"
#include "io/read_stream_vector.hpp"
#include "simple_tools.hpp"
#include <unordered_set>

namespace debruijn_graph {

extern const string tmp_folder;

//void ConstructGraphFromGenome(size_t k, Graph& g, EdgeIndex<Graph>& index/*, CoverageHandler<DeBruijnGraph>& coverage_handler*/
//, PairedInfoIndexT<Graph>& paired_index, const string& genome
//		, size_t read_size) {
//	typedef read_generator::ReadGenerator<read_generator::SmoothPositionChooser> Stream;
//	size_t coverage = 2 * read_size;
//	size_t gap = 0;
//	Stream raw_stream(2, read_size, genome, coverage, gap);
//	typedef io::RCReaderWrapper<io::PairedRead> RCStream;
//	RCStream read_stream(raw_stream);
//	ConstructGraphWithPairedInfo (k, g, index/*, coverage_handler*/,
//			paired_index, read_stream);
//}

using io::SingleRead;
using io::PairedRead;

typedef string MyRead;
typedef pair<MyRead, MyRead> MyPairedRead;
typedef string MyEdge;
typedef pair<MyEdge, MyEdge> MyEdgePair;
typedef multimap<MyEdgePair, pair<int, double>> EdgePairInfo;
typedef map<MyEdge, double> CoverageInfo;
typedef unordered_set<MyEdge> Edges;

string print(const Edges& es) {
	string s = "Edge set : {";
	for (auto i = es.begin(); i != es.end(); ++i) {
		s += "'" + *i + "'; ";
	}
	return s + "}";
}

const Edges AddComplement(const Edges& edges) {
	Edges ans;
	for (auto it = edges.begin(); it != edges.end(); ++it) {
		ans.insert(*it);
		ans.insert(ReverseComplement(*it));
	}
	return ans;
}

const CoverageInfo AddComplement(const CoverageInfo& coverage_info) {
	CoverageInfo ans;
	for (auto it = coverage_info.begin(); it != coverage_info.end(); ++it) {
		ans.insert(*it);
		ans.insert(make_pair(ReverseComplement((*it).first), (*it).second));
	}
	return ans;
}

const EdgePairInfo AddBackward(const EdgePairInfo& pair_info) {
	EdgePairInfo ans;
	for (auto it = pair_info.begin(); it != pair_info.end(); ++it) {
		ans.insert(*it);
		ans.insert(make_pair(make_pair((*it).first.second, (*it).first.first), make_pair(-(*it).second.first, (*it).second.second)));
	}
	return ans;
}

const EdgePairInfo AddComplement(const EdgePairInfo& pair_info) {
	EdgePairInfo ans;
	for (auto it = pair_info.begin(); it != pair_info.end(); ++it) {
		ans.insert(*it);
		ans.insert(make_pair(make_pair(ReverseComplement((*it).first.second), ReverseComplement((*it).first.first)), (*it).second));
	}
	return ans;
}

void EdgesEqual(const Edges& s1, const Edges& s2) {
	BOOST_CHECK_EQUAL(s2.size(), s1.size());
	for (auto it = s1.begin(); it != s1.end(); ++it) {
		BOOST_CHECK(s2.count(*it) > 0);
	}
}

const io::SingleRead MakeRead(const MyRead& read) {
	//todo fill with good quality
	std::string qual;
	qual.resize(read.size());
	return io::SingleRead("", read, qual);
}

const vector<io::SingleRead> MakeReads(const vector<MyRead>& reads) {
	vector<io::SingleRead> ans;
	for (size_t i = 0; i < reads.size(); ++i) {
		ans.push_back(MakeRead(reads[i]));
	}
	return ans;
}

const vector<PairedRead> MakePairedReads(const vector<MyPairedRead>& paired_reads, size_t insert_size) {
	DEBUG("Making paired reads");
	vector<PairedRead> ans;
	for (size_t i = 0; i < paired_reads.size(); ++i) {
		ans.push_back(PairedRead(MakeRead(paired_reads[i].first), MakeRead(paired_reads[i].second), insert_size));
	}
	DEBUG("Made paired reads");
	return ans;
}

void AssertEdges(Graph& g, const Edges& etalon_edges) {
	DEBUG("Asserting edges");
	Edges edges;
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		edges.insert(g.EdgeNucls(*it).str());
	}
	EdgesEqual(edges, etalon_edges);
}

void AssertGraph(size_t k, const vector<string>& reads, const vector<string>& etalon_edges) {
	DEBUG("Asserting graph");
	typedef io::VectorReader<SingleRead> RawStream;
	typedef io::RCReaderWrapper<SingleRead> Stream;
	RawStream raw_stream(MakeReads(reads));
	Stream read_stream(raw_stream);
	Graph g(k);
	EdgeIndex<Graph> index(g, k + 1, tmp_folder);

	io::ReadStreamVector< io::IReader<SingleRead> > streams(&read_stream);
	ConstructGraph(k, streams, g, index);

	AssertEdges(g, AddComplement(Edges(etalon_edges.begin(), etalon_edges.end())));
}

bool EqualDouble(double d1, double d2) {
	return std::abs(d1 - d2) < 1e-5;
}

void AssertCoverage(Graph& g, const CoverageInfo& etalon_coverage) {
	DEBUG("Asserting coverage");
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		auto etalon_cov_it = etalon_coverage.find(g.EdgeNucls(*it).str());
		BOOST_CHECK_MESSAGE(etalon_cov_it != etalon_coverage.end(), "Etalon didn't contain edge '" << g.EdgeNucls(*it) << "'");
		BOOST_CHECK_MESSAGE(EqualDouble(g.coverage(*it), (*etalon_cov_it).second),
				"Coverage for edge '" << g.EdgeNucls(*it) << "' was " << g.coverage(*it) << " but etalon is " << (*etalon_cov_it).second);
	}
}

typedef PairedInfoIndexT<Graph> PairedIndex;
typedef vector<PairInfo<EdgeId>> PairInfos;

void AssertPairInfo(const Graph& g, /*todo const */PairedIndex& paired_index, const EdgePairInfo& etalon_pair_info) {
	for (auto it = paired_index.begin(); it != paired_index.end(); ++it) {
		PairInfos infos;
    infos.reserve(it->size());
    for (auto set_it = it->begin(); set_it != it->end(); ++set_it)
      infos.push_back(PairInfo<EdgeId>(it.first(), it.second(), *set_it));
		for (auto info_it = infos.begin(); info_it != infos.end(); ++info_it) {
			PairInfo<EdgeId> pair_info = *info_it;
			if (pair_info.first == pair_info.second && rounded_d(pair_info) == 0) {
				continue;
			}
			pair<MyEdge, MyEdge> my_edge_pair(g.EdgeNucls(pair_info.first).str(), g.EdgeNucls(pair_info.second).str());
			auto equal_range = etalon_pair_info.equal_range(my_edge_pair);

			string my_edge_pair_str = "[" + my_edge_pair.first + ", " + my_edge_pair.second + "]";
			BOOST_CHECK_MESSAGE(equal_range.first != equal_range.second,
					"Pair of edges " << my_edge_pair_str << " wasn't found in etalon");

			double etalon_weight = -1.0;

			for (auto range_it = equal_range.first; range_it != equal_range.second; ++range_it) {
				if ((*range_it).second.first == rounded_d(pair_info)) {
					etalon_weight = (*range_it).second.second;
				}
			}
			BOOST_CHECK_MESSAGE(etalon_weight > 0,
					"Etalon didn't contain distance=" << rounded_d(pair_info) << " for edge pair " << my_edge_pair_str);
			BOOST_CHECK_MESSAGE(EqualDouble(etalon_weight, pair_info.weight()),
					"Actual weight for edge pair " << my_edge_pair_str << " on distance " << rounded_d(pair_info) << " was " << pair_info.weight() << " but etalon is " <<  etalon_weight);
		}
	}
}

void AssertGraph(size_t k, const vector<MyPairedRead>& paired_reads, size_t insert_size, const vector<MyEdge>& etalon_edges
		, const CoverageInfo& etalon_coverage, const EdgePairInfo& etalon_pair_info) {
	typedef io::VectorReader<PairedRead> RawStream;
	typedef io::RCReaderWrapper<PairedRead> PairedStream;
	typedef io::ConvertingReaderWrapper SingleStream;

	DEBUG("Asserting graph with etalon data");

	RawStream raw_stream(MakePairedReads(paired_reads, insert_size));
	PairedStream paired_read_stream(raw_stream);
	io::ReadStreamVector<io::IReader<io::PairedRead>> paired_stream_vector({&paired_read_stream});
	DEBUG("Streams initialized");

	graph_pack<Graph, runtime_k::RtSeq> gp(k, tmp_folder, (Sequence()));
	DEBUG("Graph pack created");

	PairedInfoIndexT<Graph> paired_index(gp.g);

	SingleStream single_stream(paired_read_stream);
	io::ReadStreamVector<io::IReader<io::SingleRead>> single_stream_vector({&single_stream});
	ConstructGraphWithCoverage<io::SingleRead>(k, single_stream_vector, gp.g, gp.index);

	FillPairedIndexWithReadCountMetric(gp.g,
			gp.int_ids, gp.index,
			gp.kmer_mapper,
			paired_index,
			paired_stream_vector, k);

	AssertEdges(gp.g, AddComplement(Edges(etalon_edges.begin(), etalon_edges.end())));

	AssertCoverage(gp.g, AddComplement(etalon_coverage));

	AssertPairInfo(gp.g, paired_index, AddComplement(AddBackward(etalon_pair_info)));
}

struct TmpFolderFixture
{
	TmpFolderFixture()
    {
        make_dir(tmp_folder);
    }

    ~TmpFolderFixture()
    {
    	//todo check that folder is empty
    	remove_dir(tmp_folder);
//        BOOST_TEST_MESSAGE("teardown mass");
    }
};

}
