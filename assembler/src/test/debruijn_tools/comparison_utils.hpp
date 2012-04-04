//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "graphio.hpp"
#include "xmath.h"
#include <iostream>
#include "logging.hpp"
#include "io/splitting_wrapper.hpp"
#include "io/vector_reader.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

namespace debruijn_graph {

template<size_t k>
void FillBagForStrand(const Sequence& strand, map<Seq<k>, size_t, typename Seq<k>::less2>& bag) {
	if (strand.size() < k)
		return;
	Seq<k> kmer(strand);
	kmer >> 'A';
	for (size_t i = k - 1; i < strand.size(); ++i) {
		kmer = kmer << strand[i];
		bag[kmer] += 1;
	}
}

template<size_t k>
void FillRepeats(const Sequence& genome, set<Seq<k>, typename Seq<k>::less2>& repeats) {
	map<Seq<k>, size_t, typename Seq<k>::less2> bag;

	FillBagForStrand(genome, bag);
	FillBagForStrand(!genome, bag);

	for (auto it = bag.begin(); it != bag.end(); ++it) {
		if (it->second > 1)
			repeats.insert(it->first);
	}
}

template<size_t k>
Sequence ClearGenome(const Sequence& genome, const set<Seq<k>, typename Seq<k>::less2>& repeats) {
	INFO("Clearing genome");
	if (genome.size() < k)
		return genome;

	string answer;
	for (size_t i = 0; i < k - 1; ++i) {
		answer += nucl(genome[i]);
	}
	//intervals of kmers that should be resolved afterwards
	vector<Range> repeat_intervals;
	Seq<k> kmer(genome);
	size_t curr_pos = 0;
	//curr_pos + k - next nucl pos
	bool changed = false;
	while(curr_pos + k != genome.size()) {
		size_t int_start = curr_pos;
		while (repeats.count(kmer) > 0 && curr_pos + k < genome.size()) {
			kmer = kmer << genome[curr_pos + k];
			curr_pos++;
			changed = true;
		}

		repeat_intervals.push_back(Range(int_start, curr_pos));

		if (curr_pos + k == genome.size())
			break;

		while (repeats.count(kmer) == 0 && curr_pos + k < genome.size()) {
			answer += nucl(kmer[k - 1]);
			kmer = kmer << genome[curr_pos + k];
			curr_pos++;
		}
	}
	if (changed) {
		INFO("Genome was changed during cleaning");
	} else {
		INFO("Genome wasn't changed during cleaning");
	}
	return Sequence(answer);
}

template<size_t k>
Sequence ClearGenome(const Sequence& genome) {
	INFO("Clearing genome of repeats");

	set<Seq<k>, typename Seq<k>::less2> repeats;
	INFO("Filling set of repeats");
	FillRepeats<k>(genome, repeats);
	INFO("Clearing genome");
	return ClearGenome<k>(genome, repeats);
}

//todo bad strategy for assembly cleaning
template<size_t k>
pair<Sequence, vector<Sequence>> Clear(const Sequence& genome, const vector<Sequence>& assembly) {
	INFO("Clearing genome of repeats");

	set<Seq<k>, typename Seq<k>::less2> repeats;
	INFO("Filling set of repeats");
	FillRepeats<k>(genome, repeats);
	for (auto it = assembly.begin(); it != assembly.end(); ++it) {
		FillRepeats(*it, repeats);
	}
	INFO("Clearing genome");
	Sequence new_genome = ClearGenome<k>(genome, repeats);
	INFO("Clearing assembly");
	vector<Sequence> new_assembly;
	for (auto it = assembly.begin(); it != assembly.end(); ++it) {
		new_assembly.push_back(ClearGenome<k>(*it, repeats));
	}
	return make_pair(new_genome, new_assembly);
}

template<size_t k>
pair<Sequence, Sequence> ClearGenomes(const pair<Sequence, Sequence>& genomes) {
	INFO("Clearing genomes from repeats");

	set<Seq<k>, typename Seq<k>::less2> repeats;
	INFO("Filling set of repeats");
	FillRepeats<k>(genomes.first, repeats);
	FillRepeats<k>(genomes.second, repeats);
	INFO("Clearing genomes");
	return make_pair(ClearGenome<k>(genomes.first, repeats), ClearGenome<k>(genomes.second, repeats));
}

template<size_t k>
pair<Sequence, Sequence> TotallyClearGenomes(const pair<Sequence, Sequence>& genomes) {
	static const size_t iter_count = 1;
	pair<Sequence, Sequence> tmp = genomes;
	for (size_t i = 0; i < iter_count; ++i) {
		INFO("Cleaning iteration " << i);
		tmp = ClearGenomes<k>(tmp);
	}
	return tmp;
}

double uniform_01() {
	static boost::mt19937 rng(43);
	static boost::uniform_01<boost::mt19937> zeroone(rng);
	return zeroone();
}

bool event_happened(double rate) {
	return ls(uniform_01(), rate);
}

int rand_int(size_t min, size_t max) {
	static boost::mt19937 rng(43);
	boost::uniform_int<> un_int(min, max);
	boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
	         die(rng, un_int);
	return die();
}

char switch_nucl(char n) {
	VERIFY(is_nucl(n));
	return nucl((dignucl(n) + rand_int(1,3)) % 4);
}

Sequence IntroduceReversal(const Sequence& s, size_t min_len, size_t max_len) {
	VERIFY(s.size() > min_len);
	//inclusive
	size_t start = rand_int(0, s.size() - min_len);
	size_t len = rand_int(min_len, std::min(max_len, s.size() - start));
	//exclusive
	size_t end = start + len;
	INFO("Reversing fragment of length " << len << " from " << start << " to " << end);
	return s.Subseq(0, start) + !s.Subseq(start, end) + s.Subseq(end);
}

Sequence IntroduceReversals(const Sequence& s, size_t rev_count, size_t min_len, size_t max_len) {
	Sequence res = s;
	for (size_t i = 0; i < rev_count; ++i) {
		res = IntroduceReversal(res, min_len, max_len);
	}
	return res;
}

template<class gp_t>
void ConstructRepeatGraph(gp_t& gp) {
	io::VectorReader<io::SingleRead> stream(
			io::SingleRead("genome", gp.genome.str()));
	io::RCReaderWrapper<io::SingleRead> rc_stream(stream);
	ConstructGraph<gp_t::k_value, typename gp_t::graph_t>(gp.g, gp.index, rc_stream);
}

template<class Graph>
vector<Sequence> EdgesSequences(const Graph& g) {
	vector<Sequence> res;
	set<EdgeId> edges;
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if (edges.find(*it) == edges.end()) {
			res.push_back(g.EdgeNucls(*it));
			edges.insert(g.conjugate(*it));
		}
	}
	return res;
}

template<class gp_t>
vector<Sequence> RepeatGraphEdges(const Sequence& genome) {
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;

	gp_t gp(genome);
	ConstructRepeatGraph(gp);
	return EdgesSequences(gp.g);
}

Sequence ReadGenome(const string& filename) {
	checkFileExistenceFATAL(filename);
	io::Reader genome_stream(filename);
	io::SingleRead genome;
	genome_stream >> genome;
	return genome.sequence();
}

Sequence IntroduceMutations(const Sequence& s, double rate) {
	VERIFY(ge(rate, 0.) && ls(rate, 1.0));
	string as_str = s.str();
	for (size_t i = 0; i < s.size(); ++i) {
		if (event_happened(rate)) {
			as_str[i] = switch_nucl(as_str[i]);
		}
	}
	return Sequence(as_str);
}

template<size_t k, size_t K>
void RunBPComparison(ContigStream& raw_stream1,
		ContigStream& raw_stream2, const string& name1, const string& name2,
		bool refine, bool untangle, const string& output_folder,
		bool detailed_output = true, size_t delta = 5) {
	io::SplittingWrapper stream1(raw_stream1);
	io::SplittingWrapper stream2(raw_stream2);

	typedef graph_pack</*Nonc*/ConjugateDeBruijnGraph, K> comparing_gp_t;
	if (refine) {
		typedef graph_pack<ConjugateDeBruijnGraph, k> refining_gp_t;
		refining_gp_t refining_gp;
		ConstructGPForRefinement(refining_gp, stream1, stream2, delta);

		ContigRefiner<refining_gp_t> refined_stream1(stream1, refining_gp);
		ContigRefiner<refining_gp_t> refined_stream2(stream2, refining_gp);

		AssemblyComparer<comparing_gp_t> comparer(refined_stream1,
				refined_stream2, name1, name2, untangle);
		comparer.CompareAssemblies(output_folder, detailed_output, /*one_many_resolve*/false);
	} else {
		AssemblyComparer<comparing_gp_t> comparer(stream1,
				stream2, name1, name2, untangle);
		comparer.CompareAssemblies(output_folder, detailed_output, /*one_many_resolve*/false);
	}
}

template<size_t k, size_t K>
void RunBPComparison(const Sequence& ref,
		ContigStream& stream, const string& name1, const string& name2,
		bool refine, bool untangle, const string& output_folder,
		bool detailed_output = true, size_t delta = 5) {
	io::VectorReader<io::SingleRead> ref_stream(
			io::SingleRead(name1, ref.str()));
	RunBPComparison<k, K>(ref_stream, stream, name1, name2,
			refine, untangle, output_folder,
			detailed_output, delta);
}

template<size_t k, size_t K>
void RunBPComparison(const Sequence& s1,
		const Sequence& s2, const string& name1, const string& name2,
		bool refine, bool untangle, const string& output_folder,
		bool detailed_output = true) {
	io::VectorReader<io::SingleRead> stream(
			io::SingleRead(name2, s2.str()));
	RunBPComparison<k, K>(s1, stream, name1, name2,
			refine, untangle, output_folder,
			detailed_output);
}

const vector<io::SingleRead> MakeReads(const vector<Sequence>& ss) {
	vector<io::SingleRead> ans;
	for (size_t i = 0; i < ss.size(); ++i) {
		ans.push_back(io::SingleRead("read_" + ToString(i), ss[i].str()));
	}
	return ans;
}

template<size_t k, size_t K>
void RunBPComparison(const Sequence& ref,
		const vector<Sequence>& contigs, const string& name1, const string& name2,
		bool refine, bool untangle, const string& output_folder,
		bool detailed_output = true) {
	io::VectorReader<io::SingleRead> stream(MakeReads(contigs));
	RunBPComparison<k, K>(ref, stream, name1, name2,
			refine, untangle, output_folder,
			detailed_output);
}

Sequence FirstSequence(io::IReader<io::SingleRead>& stream) {
	stream.reset();
	io::SingleRead r;
	VERIFY(!stream.eof());
	stream >> r;
	return r.sequence();
}

vector<Sequence> AllSequences(io::IReader<io::SingleRead>& stream) {
	vector<Sequence> answer;
	stream.reset();
	io::SingleRead r;
	while (!stream.eof()) {
		stream >> r;
		answer.push_back(r.sequence());
	}
	return answer;
}

template<size_t k>
pair<Sequence, Sequence> CorrectGenomes(const Sequence& genome1, const Sequence& genome2, size_t delta = 5) {
	io::VectorReader<io::SingleRead> stream1(
			io::SingleRead("first", genome1.str()));
	io::VectorReader<io::SingleRead> stream2(
			io::SingleRead("second", genome2.str()));

	typedef graph_pack<ConjugateDeBruijnGraph, k> refining_gp_t;
	refining_gp_t refining_gp;
	ConstructGPForRefinement(refining_gp, stream1, stream2, delta);

	ContigRefiner<refining_gp_t> refined_stream1(stream1, refining_gp);
	ContigRefiner<refining_gp_t> refined_stream2(stream2, refining_gp);

	pair<Sequence, Sequence> answer = make_pair(FirstSequence(refined_stream1), FirstSequence(refined_stream2));
	return answer;
}

template<size_t k>
pair<Sequence, Sequence> CorrectGenomes(const pair<Sequence, Sequence>& genomes, size_t delta = 5) {
	return CorrectGenomes<k>(genomes.first, genomes.second, delta);
}

template<size_t k>
bool CheckNoRepeats(const Sequence& genome) {
	set<Seq<k>, typename Seq<k>::less2> repeats;
	FillRepeats<k>(genome, repeats);
	return repeats.empty();
}

//Prints only basic graph structure!!!
//todo rewrite with normal splitter usage instead of filtering
void PrintGraphComponentContainingEdge(const string& file_name,
		const Graph& g, size_t split_edge_length,
		const IdTrackHandler<Graph>& int_ids, int int_edge_id) {
	LongEdgesInclusiveSplitter<Graph> inner_splitter(g, split_edge_length);

//	VERIFY_MSG(int_ids.ReturnEdgeId(int_edge_id) != NULL,
//			"Couldn't find edge with id = " << int_edge_id);

	AnyEdgeContainFilter<Graph> filter(g, int_ids.ReturnEdgeId(int_edge_id));
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, filter);
	vector<vector<VertexId>> components;
	while (!splitter.Finished()) {
		components.push_back(splitter.NextComponent());
	}VERIFY(components.size() == 1);
	ConjugateDataPrinter<Graph> printer(g, components.front().begin(),
			components.front().end(), int_ids);
	PrintBasicGraph<Graph>(file_name, printer);
}

template<class gp_t>
void ThreadAssemblies(const string& base_saves,
		ContigStream& base_assembly, const string& base_prefix,
		ContigStream& assembly_to_thread, const string& to_thread_prefix,
		const string& output_dir) {
	typedef typename gp_t::graph_t Graph;
	gp_t gp;
//		ConstructGraph<gp_t::k_value, Graph>(gp.g, gp.index, base_assembly);
	ScanGraphPack(base_saves, gp);
	base_assembly.reset();
	FillPos(gp, base_assembly, base_prefix);
	FillPos(gp, assembly_to_thread, to_thread_prefix);

	EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);
	StrGraphLabeler<Graph> str_labeler(gp.g);
	CompositeLabeler<Graph> labeler(pos_labeler, str_labeler);

	NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> mapper(gp.g, gp.index,
			gp.kmer_mapper);

	assembly_to_thread.reset();
	io::SingleRead read;
	while (!assembly_to_thread.eof()) {
		assembly_to_thread >> read;
		make_dir(output_dir + read.name());
		WriteComponentsAlongPath(gp.g, labeler,
				output_dir + read.name() + "/.dot", /*split_edge_length*/400,
				mapper.MapSequence(read.sequence()),
				Path<typename Graph::EdgeId>(), Path<typename Graph::EdgeId>(),
				true);
	}
}

template <size_t k>
pair<Sequence, vector<Sequence>> RefineData(const pair<Sequence, vector<Sequence>>& data) {
	io::VectorReader<io::SingleRead> stream1(
			io::SingleRead("first", data.first.str()));
	io::VectorReader<io::SingleRead> stream2(MakeReads(data.second));

	typedef graph_pack<ConjugateDeBruijnGraph, k> refining_gp_t;
	refining_gp_t refining_gp;
	ConstructGPForRefinement(refining_gp, stream1, stream2);

	ContigRefiner<refining_gp_t> refined_stream1(stream1, refining_gp);
	ContigRefiner<refining_gp_t> refined_stream2(stream2, refining_gp);

	return make_pair(FirstSequence(refined_stream1), AllSequences(refined_stream2));
}

template <size_t k>
void CompareGenomes(const Sequence& genome_1, const Sequence& genome_2, const string& output_dir) {
	INFO("Genome comparison started");
	io::VectorReader<io::SingleRead> stream1(
			io::SingleRead("", genome_1.str()));
	io::VectorReader<io::SingleRead> stream2(
			io::SingleRead("", genome_2.str()));
	typedef graph_pack</*Nonc*/ConjugateDeBruijnGraph, k> comparing_gp_t;
	INFO("Running assembly comparer");
	AssemblyComparer<comparing_gp_t> comparer(stream1, stream2, "strain1", "strain2", /*untangle*/false);
	comparer.CompareAssemblies(output_dir, /*detailed_output*/true, /*on_many_resolve*/true);
	INFO("Finished");
}

}
