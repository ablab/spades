//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "graphio.hpp"
#include "simple_tools.hpp"
#include "new_debruijn.hpp"
#include "xmath.h"
#include <iostream>
#include "logger/logger.hpp"
#include "io/multifile_reader.hpp"
#include "io/splitting_wrapper.hpp"
#include "io/modifying_reader_wrapper.hpp"
#include "io/vector_reader.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

namespace compare {
using namespace debruijn_graph;

inline Sequence ReadSequence(io::IReader<io::SingleRead>& reader) {
	VERIFY(!reader.eof());
	io::SingleRead read;
	reader >> read;
	return read.sequence();
}

template<size_t k, class Graph>
void ConstructGraph(Graph& g, EdgeIndex<k + 1, Graph>& index,
		io::IReader<io::SingleRead>& stream) {
	vector<io::IReader<io::SingleRead>*> streams = { &stream };
	ConstructGraph<k, Graph>(streams, g, index);
}

template<size_t k, class Graph>
void ConstructGraph(Graph& g, EdgeIndex<k + 1, Graph>& index,
		io::IReader<io::SingleRead>& stream1,
		io::IReader<io::SingleRead>& stream2) {
	io::MultifileReader<io::SingleRead> composite_reader(stream1, stream2);
	ConstructGraph<k, Graph>(g, index, composite_reader);
}

inline Sequence ReadGenome(const string& filename) {
	checkFileExistenceFATAL(filename);
	io::Reader genome_stream(filename);
	return ReadSequence(genome_stream);
}

inline vector<io::SingleRead> MakeReads(const vector<Sequence>& ss) {
	vector<io::SingleRead> ans;
	for (size_t i = 0; i < ss.size(); ++i) {
		ans.push_back(io::SingleRead("read_" + ToString(i), ss[i].str()));
	}
	return ans;
}

inline Sequence FirstSequence(io::IReader<io::SingleRead>& stream) {
	stream.reset();
	io::SingleRead r;
	VERIFY(!stream.eof());
	stream >> r;
	return r.sequence();
}

inline vector<Sequence> AllSequences(io::IReader<io::SingleRead>& stream) {
	vector<Sequence> answer;
	stream.reset();
	io::SingleRead r;
	while (!stream.eof()) {
		stream >> r;
		answer.push_back(r.sequence());
	}
	return answer;
}

inline vector<Sequence> ReadContigs(const string& filename) {
	checkFileExistenceFATAL(filename);
	io::Reader genome_stream(filename);
	return AllSequences(genome_stream);
}

//Prints only basic graph structure!!!
//todo rewrite with normal splitter usage instead of filtering
inline void PrintGraphComponentContainingEdge(const string& file_name, const Graph& g,
		size_t split_edge_length, const IdTrackHandler<Graph>& int_ids,
		int int_edge_id) {
	LongEdgesInclusiveSplitter<Graph> inner_splitter(g, split_edge_length);

//	VERIFY_MSG(int_ids.ReturnEdgeId(int_edge_id) != NULL,
//			"Couldn't find edge with id = " << int_edge_id);

	AnyEdgeContainFilter<Graph> filter(g, int_ids.ReturnEdgeId(int_edge_id));
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, filter);
	vector<vector<VertexId>> components;
	while (!splitter.Finished()) {
		components.push_back(splitter.NextComponent());
	}
	VERIFY(components.size() == 1);
	ConjugateDataPrinter<Graph> printer(g, components.front().begin(),
			components.front().end(), int_ids);
	PrintBasicGraph<Graph>(file_name, printer);
}

}
