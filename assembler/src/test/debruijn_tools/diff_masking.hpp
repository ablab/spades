#pragma once

#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "graph_read_correction.hpp"

namespace compare {

template<class gp_t>
void ConstructGPForRefinement(gp_t& gp, const vector<ContigStream*>& contigs,
		size_t delta) {
	typedef typename gp_t::graph_t Graph;
	INFO("Constructing graph pack for refinement");

	vector<ContigStream*> rc_contigs;
	for (auto it = contigs.begin(); it != contigs.end(); ++it) {
		rc_contigs.push_back(new RCWrapper(**it));
		rc_contigs.back()->reset();
	}

	CompositeContigStream comp_stream(rc_contigs);

	ConstructGraph<gp_t::k_value>(gp.g, gp.index, comp_stream);

//	make_dir("bp_graph_test/tmp/");
//	LengthIdGraphLabeler<Graph> labeler(gp.g);
//	WriteToDotFile(gp.g, labeler, "bp_graph_test/tmp/before_refine.dot");

//todo configure!!!
	debruijn_config::simplification::bulge_remover br_config;
	br_config.max_bulge_length_coefficient = 3;
	br_config.max_coverage = 1000.;
	br_config.max_relative_coverage = 1.2;
	br_config.max_delta = delta;
	br_config.max_relative_delta = 0.1;

	INFO("Removing bulges");
	RemoveBulges(gp.g, br_config);

	INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");

	TipsProjector<gp_t> tip_projector(gp);
	boost::function<void(EdgeId)> projecting_callback = boost::bind(
			&TipsProjector<gp_t>::ProjectTip, &tip_projector, _1);
	debruijn_config::simplification::tip_clipper tc_config;
	tc_config.max_coverage = 1000.;
	tc_config.max_relative_coverage = 1.1;
	tc_config.max_tip_length_coefficient = 2.;

	INFO("Clipping tips with projection");
	DefaultClipTips(gp.g, tc_config, /*read_length*/100, projecting_callback);

	INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");

//	WriteToDotFile(gp.g, labeler, "bp_graph_test/tmp/after_refine.dot");
}

template<class gp_t>
void ConstructGPForRefinement(gp_t& gp,
		io::IReader<io::SingleRead>& raw_stream_1,
		io::IReader<io::SingleRead>& raw_stream_2, size_t delta) {
	vector<ContigStream*> contigs = { &raw_stream_1, &raw_stream_2 };
	ConstructGPForRefinement(gp, contigs, delta);
}


template<size_t k>
pair<Sequence, Sequence> CorrectGenomes(const Sequence& genome1,
		const Sequence& genome2, size_t delta = 5) {
	io::VectorReader<io::SingleRead> stream1(
			io::SingleRead("first", genome1.str()));
	io::VectorReader<io::SingleRead> stream2(
			io::SingleRead("second", genome2.str()));

	typedef graph_pack<ConjugateDeBruijnGraph, k> refining_gp_t;
	refining_gp_t refining_gp;
	ConstructGPForRefinement(refining_gp, stream1, stream2, delta);

	io::ModifyingWrapper<io::SingleRead> refined_stream1(stream1
			, GraphReadCorrectorInstance(refining_gp.g, *MapperInstance(refining_gp)));
	io::ModifyingWrapper<io::SingleRead> refined_stream2(stream2
			, GraphReadCorrectorInstance(refining_gp.g, *MapperInstance(refining_gp)));

	pair<Sequence, Sequence> answer = make_pair(FirstSequence(refined_stream1),
			FirstSequence(refined_stream2));
	return answer;
}

template<size_t k>
pair<Sequence, Sequence> CorrectGenomes(const pair<Sequence, Sequence>& genomes
		, size_t delta = 5) {
	return CorrectGenomes<k>(genomes.first, genomes.second, delta);
}


template<size_t k>
pair<Sequence, vector<Sequence>> RefineData(
		const pair<Sequence, vector<Sequence>>& data) {
	io::VectorReader<io::SingleRead> stream1(
			io::SingleRead("first", data.first.str()));
	io::VectorReader<io::SingleRead> stream2(MakeReads(data.second));

	typedef graph_pack<ConjugateDeBruijnGraph, k> refining_gp_t;
	refining_gp_t refining_gp;
	ConstructGPForRefinement(refining_gp, stream1, stream2);

	io::ModifyingWrapper<io::SingleRead> refined_stream1(stream1
			, GraphReadCorrectorInstance(refining_gp.g, *MapperInstance(refining_gp)));
	io::ModifyingWrapper<io::SingleRead> refined_stream2(stream2
			, GraphReadCorrectorInstance(refining_gp.g, *MapperInstance(refining_gp)));

	return make_pair(FirstSequence(refined_stream1),
			AllSequences(refined_stream2));
}

}
