//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "io/read_stream_vector.hpp"
#include "graph_construction.hpp"
#include "graph_simplification.hpp"
#include "graph_read_correction.hpp"
#include "test_utils.hpp"

#include "coloring.hpp"
#include "colored_graph_construction.hpp"

namespace cap {

template<class Stream1, class Stream2>
void Transfer(Stream1& s1, Stream2& s2) {
    typename Stream1::read_type r;
    while (!s1.eof()) {
        s1 >> r;
        s2 << r;
    }
}

//todo changes the graph!!! color edge splitting!!!
template<class gp_t>
void MakeSaves(gp_t& gp, ContigStreamsPtr streams, const string& root,
               const vector<string>& out_files, bool optional = true) {

    static bool make_optional_saves = true;

    if (!make_optional_saves && optional)
        return;

    make_dir(root);

    using namespace debruijn_graph;
    ContigStreamsPtr rc_contigs = io::RCWrapStreams(*streams);

    rc_contigs->reset();

    ColorHandler<Graph> coloring(gp.g, rc_contigs->size());
    SplitAndColorGraph(gp, coloring, *rc_contigs, true);

    PrintColoredGraphWithColorFilter(gp.g, coloring, gp.edge_pos,
                                     root + "colored_split_graph.dot");
}

template<class gp_t>
void RefineGP(gp_t& gp, size_t delta = 5) {
    using namespace debruijn_graph;
    typedef typename gp_t::graph_t Graph;
    INFO("Constructing graph pack for refinement");

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

    debruijn_config::simplification::complex_bulge_remover cbr_config;
    cbr_config.enabled = true;
    cbr_config.pics_enabled = false;
    cbr_config.folder = "";
    cbr_config.max_relative_length = 3;
    cbr_config.max_length_difference = 1000;

    INFO("Removing complex bulges");
    RemoveComplexBulges(gp.g, cbr_config);

    TipsProjector<gp_t> tip_projector(gp);
    boost::function<void(EdgeId)> projecting_callback = boost::bind(
            &TipsProjector<gp_t>::ProjectTip, &tip_projector, _1);
    debruijn_config::simplification::tip_clipper tc_config;

    tc_config.condition = "{ tc_lb 2. }";

    INFO("Clipping tips with projection");

    ClipTipsWithProjection(gp, tc_config, true);

    INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");
}

template<class gp_t>
void ConstructGPForRefinement(gp_t& gp, const ContigStreamsPtr& contigs,
                              size_t delta = 5) {
    using namespace debruijn_graph;
    typedef typename gp_t::graph_t Graph;
    INFO("Constructing graph pack for refinement");

    ContigStreamsPtr rc_streams = io::RCWrapStreams(*contigs);
    rc_streams->reset();

    ConstructGraph(gp.k_value, *rc_streams, gp.g, gp.index);

    RefineGP(gp, delta);
}

template<class gp_t>
void ConstructGPForRefinement(gp_t& gp,
                              io::IReader<io::SingleRead>& raw_stream_1,
                              io::IReader<io::SingleRead>& raw_stream_2,
                              size_t delta = 5) {
    ContigStreamsPtr streams_ptr = make_shared<ContigStreams>(
        vector<ContigStream*>{&raw_stream_1, &raw_stream_2}, false);
    ConstructGPForRefinement(gp, streams_ptr, delta);
}

template<size_t k, class Seq>
pair<Sequence, Sequence> CorrectGenomes(const Sequence& genome1,
                                        const Sequence& genome2, size_t delta =
                                                5) {
    io::VectorReader<io::SingleRead> stream1(
            io::SingleRead("first", genome1.str()));
    io::VectorReader<io::SingleRead> stream2(
            io::SingleRead("second", genome2.str()));

    typedef debruijn_graph::graph_pack<debruijn_graph::ConjugateDeBruijnGraph,
            Seq> refining_gp_t;
    refining_gp_t refining_gp(k, "tmp");
    ConstructGPForRefinement(refining_gp, stream1, stream2, delta);

    io::ModifyingWrapper<io::SingleRead> refined_stream1(
            stream1,
            GraphReadCorrectorInstance(refining_gp.g,
                                       *MapperInstance(refining_gp)));
    io::ModifyingWrapper<io::SingleRead> refined_stream2(
            stream2,
            GraphReadCorrectorInstance(refining_gp.g,
                                       *MapperInstance(refining_gp)));

    pair<Sequence, Sequence> answer = make_pair(FirstSequence(refined_stream1),
                                                FirstSequence(refined_stream2));
    return answer;
}

template<size_t k>
pair<Sequence, Sequence> CorrectGenomes(const pair<Sequence, Sequence>& genomes,
                                        size_t delta = 5) {
    return CorrectGenomes<k>(genomes.first, genomes.second, delta);
}

template<size_t k, class Seq>
pair<Sequence, vector<Sequence>> RefineData(
        const pair<Sequence, vector<Sequence>>& data) {
    io::VectorReader<io::SingleRead> stream1(
            io::SingleRead("first", data.first.str()));
    io::VectorReader<io::SingleRead> stream2(MakeReads(data.second));

    typedef graph_pack<ConjugateDeBruijnGraph, Seq> refining_gp_t;
    refining_gp_t refining_gp(k, "tmp");
    ConstructGPForRefinement(refining_gp, stream1, stream2);

    io::ModifyingWrapper<io::SingleRead> refined_stream1(
            stream1,
            GraphReadCorrectorInstance(refining_gp.g,
                                       *MapperInstance(refining_gp)));
    io::ModifyingWrapper<io::SingleRead> refined_stream2(
            stream2,
            GraphReadCorrectorInstance(refining_gp.g,
                                       *MapperInstance(refining_gp)));

    return make_pair(FirstSequence(refined_stream1),
                     AllSequences(refined_stream2));
}

inline ContigStreamsPtr OpenStreams(const string& root, const vector<string>& filenames) {
    ContigStreamsPtr streams(new ContigStreams());
    FOREACH (auto filename, filenames) {
        DEBUG("Opening stream from " << root << filename);
        streams->push_back(new io::Reader(root + filename));
    }
    return streams;
}

inline void SaveAll(ContigStreamsPtr streams,
                    const vector<string>& suffixes, const string& out_root) {
    make_dir(out_root);

    streams->reset();
    for (size_t i = 0; i < streams->size(); ++i) {
        if (!suffixes[i].empty()) {
            string output_filename = out_root + suffixes[i] + ".fasta";
            io::ofastastream out_stream(output_filename);
            Transfer((*streams)[i], out_stream);
        }
    }
}

template<class Seq>
void MaskDifferencesAndSave(ContigStreamsPtr streams, const string& root,
                            const vector<string>& out_files, size_t k) {
    VERIFY(streams->size() == out_files.size());

    const size_t delta = std::max(size_t(5), k / 5);
    typedef graph_pack<ConjugateDeBruijnGraph, Seq> gp_t;
    typedef NewExtendedSequenceMapper<Graph, typename gp_t::seq_t> Mapper;

    make_dir(root);
    INFO("Constructing graph pack for k=" << k << " delta=" << delta);
    gp_t gp(k, "tmp");

    ContigStreamsPtr rc_streams = io::RCWrapStreams(*streams);
    rc_streams->reset();

    ConstructGraph(gp.k_value, *rc_streams, gp.g, gp.index);

    MakeSaves(gp, streams, root + "before_refinement/", out_files);

    RefineGP(gp, delta);

    ContigStreamsPtr refined_streams(new ContigStreams());
    for (size_t i = 0; i < streams->size(); ++i) {
        string output_filename = out_files[i];
        if (!output_filename.empty()) {
            refined_streams->push_back(
                    new io::ModifyingWrapper<io::SingleRead>(
                            (*streams)[i],
                            GraphReadCorrectorInstance(gp.g,
                                                       *MapperInstance(gp))));
        }
    }

    MakeSaves(gp, refined_streams, root + "after_refinement/", out_files);
}

inline void MaskDifferencesAndSave(ContigStreamsPtr streams,
                                   const vector<string>& suffixes,
                                   const string& out_root) {
    SaveAll(streams, suffixes, out_root + "final_contigs");
}

inline void MaskDifferencesAndSave(ContigStreamsPtr streams,
                                   const vector<string>& suffixes,
                                   const string& out_root,
                                   vector<size_t> &k_values) {

    if (k_values.size() == 0) {
        MaskDifferencesAndSave(streams, suffixes, out_root);
        return;
    }

    size_t current_k = k_values.back();
    k_values.pop_back();

    string root = out_root + ToString(current_k);

    if (utils::NeedToUseLongSeq(current_k)) {
        omp_set_num_threads(1);
        MaskDifferencesAndSave<LSeq>(
                streams, out_root, suffixes,
                current_k);
    } else {
        omp_set_num_threads(8);
        MaskDifferencesAndSave<runtime_k::RtSeq>(
                streams, out_root, suffixes,
                current_k);
    }

    ContigStreamsPtr corr_streams = OpenStreams(out_root, suffixes);
    //recursive call
    MaskDifferencesAndSave(corr_streams, suffixes, out_root, k_values);

}

inline void MaskDifferencesAndSave(const vector<string>& in_files,
                                   const vector<string>& suffixes,
                                   const string& out_root,
                                   vector<size_t>& k_values) {
//	remove_dir(out_root);
    utils::MakeDirPath(out_root);
    ContigStreamsPtr streams = OpenStreams("", in_files);
    MaskDifferencesAndSave(streams, suffixes, out_root, k_values);
}

}
