#pragma once

#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "graph_read_correction.hpp"

namespace cap {

template<class Collection>
void DisposeCollection(Collection c) {
	for (auto it = c.begin(); it != c.end(); ++it) {
		delete *it;
	}
}

template<class Stream1, class Stream2>
void Transfer(Stream1& s1, Stream2& s2) {
	typename Stream1::read_type r;
	while (!s1.eof()) {
		s1 >> r;
		s2 << r;
	}
}

template<class gp_t>
void ConstructGPForRefinement(gp_t& gp, const vector<ContigStream*>& contigs,
		size_t delta = 5) {
	using namespace debruijn_graph;
	typedef typename gp_t::graph_t Graph;
	INFO("Constructing graph pack for refinement");

	vector<ContigStream*> rc_contigs;
	for (auto it = contigs.begin(); it != contigs.end(); ++it) {
		rc_contigs.push_back(new RCWrapper(**it));
		rc_contigs.back()->reset();
	}

	io::ReadStreamVector<ContigStream> rc_read_stream_vector(rc_contigs);

	ConstructGraph(gp.k_value, rc_read_stream_vector, gp.g, gp.index);

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

	auto tip_clipper_factory = GetDefaultTipClipperFactory<Graph>(tc_config, /*read_length*/10000, projecting_callback);
	ClipTips(gp.g, tip_clipper_factory);

	INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");

	for (auto it = rc_contigs.begin(); it != rc_contigs.end(); ++it) {
		delete *it;
	}
//	WriteToDotFile(gp.g, labeler, "bp_graph_test/tmp/after_refine.dot");
}

template<class gp_t>
void ConstructGPForRefinement(gp_t& gp,
		io::IReader<io::SingleRead>& raw_stream_1,
		io::IReader<io::SingleRead>& raw_stream_2, size_t delta = 5) {
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

	typedef debruijn_graph::graph_pack<debruijn_graph::ConjugateDeBruijnGraph> refining_gp_t;
	refining_gp_t refining_gp(k);
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

	typedef graph_pack<ConjugateDeBruijnGraph> refining_gp_t;
	refining_gp_t refining_gp(k);
	ConstructGPForRefinement(refining_gp, stream1, stream2);

	io::ModifyingWrapper<io::SingleRead> refined_stream1(stream1
			, GraphReadCorrectorInstance(refining_gp.g, *MapperInstance(refining_gp)));
	io::ModifyingWrapper<io::SingleRead> refined_stream2(stream2
			, GraphReadCorrectorInstance(refining_gp.g, *MapperInstance(refining_gp)));

	return make_pair(FirstSequence(refined_stream1),
			AllSequences(refined_stream2));
}

template<size_t k>
void MaskDifferencesAndSave(/*const */vector<ContigStream*>& streams, const vector<string> out_files) {
	VERIFY(streams.size() == out_files.size());
	const size_t delta = std::max(size_t(5), k / 5);
	typedef graph_pack<ConjugateDeBruijnGraph> gp_t;
	DEBUG("Constructing graph pack for k=" << k << " delta=" << delta);
	gp_t gp(k);
	ConstructGPForRefinement(gp, streams, delta);
	for (size_t i = 0; i < streams.size(); ++i) {
		string output_filename = out_files[i];
		if (!output_filename.empty()) {
			io::ModifyingWrapper<io::SingleRead> refined_stream(*streams[i]
					, GraphReadCorrectorInstance(gp.g, *MapperInstance(gp)));
			refined_stream.reset();
			Contig contig;
			io::ofastastream out_stream(output_filename);
			DEBUG("Saving to " << output_filename);
			while (!refined_stream.eof()) {
				refined_stream >> contig;
				out_stream << contig;
			}
		}
	}
}

inline vector<string> CorrectPaths(const vector<string>& out_files_suffs, const string& out_root, size_t k) {
	vector<string> answer;
	for (auto it = out_files_suffs.begin(); it != out_files_suffs.end(); ++it) {
		if (it->empty()) {
			answer.push_back("");
		} else {
			stringstream ss;
			ss << out_root << k << "/" << *it << ".fasta";
			answer.push_back(ss.str());
		}
	}
	return answer;
}

inline vector<ContigStream*> OpenStreams(const vector<string>& filenames) {
	vector<ContigStream*> streams;
	for (auto it = filenames.begin(); it != filenames.end(); ++it) {
		DEBUG("Opening stream from " << *it);
		streams.push_back(new io::Reader(*it));
	}
	return streams;
}

inline void MaskDifferencesAndSave(/*const */vector<ContigStream*>& streams, const vector<string>& suffixes, const string& out_root) {
	for (size_t i = 0; i < streams.size(); ++i) {
		if (!suffixes[i].empty()) {
			string output_filename = out_root + suffixes[i] + ".fasta";
			io::ofastastream out_stream(output_filename);
			Transfer(*streams[i], out_stream);
		}
	}
}

template<size_t k_value>
struct k
{
   enum { value = k_value };
};
//
//void foo()
//{
//   cout << "final foo" << endl;
//}
//
//template <typename k, typename... ks>
//void foo(k current_k, ks... other_k)
//{
//   foo(other_k...);
//}
//
//
//
//int main()
//{
//
//  foo(k<1>(), k<2>(), k<3>());
//
//  return 0;
//}

template<typename k, typename ... Ks>
void MaskDifferencesAndSave(/*const */vector<ContigStream*>& streams, const vector<string>& suffixes
		, const string& out_root, k current_k, Ks... other_k) {
	make_dir(out_root + ToString((size_t)k::value));
	MaskDifferencesAndSave<k::value>(streams, CorrectPaths(suffixes, out_root, k::value));
	vector<ContigStream*> corr_streams = OpenStreams(CorrectPaths(suffixes, out_root, k::value));
	MaskDifferencesAndSave(corr_streams, suffixes, out_root, other_k...);

	for (auto it = corr_streams.begin(); it != corr_streams.end(); ++it) {
		delete *it;
	}
}

template<typename ... Ks>
void MaskDifferencesAndSave(const vector<string>& in_files, const vector<string>& suffixes
		, const string& out_root, Ks... ks) {
	remove_dir(out_root);
	make_dir(out_root);
	vector<ContigStream*> streams = OpenStreams(in_files);
	MaskDifferencesAndSave(streams, suffixes, out_root, ks...);
	DisposeCollection(streams);
}

//template<size_t k>
//void MaskDifferencesAndSave(/*const */vector<ContigStream*>& streams, const vector<string>& suffixes, const string& out_root) {
//	make_dir(out_root + ToString(k));
//	MaskDifferencesAndSave<k>(streams, CorrectPaths(suffixes, out_root, k));
//}
//
//template<size_t k>
//void MaskDifferencesAndSave(const vector<string>& in_files, const vector<string>& suffixes
//		, const string& out_root) {
//	make_dir(out_root);
//	vector<ContigStream*> streams = OpenStreams(in_files);
//	MaskDifferencesAndSave<k>(streams, suffixes, out_root);
//	DisposeCollection(streams);
//}

}
