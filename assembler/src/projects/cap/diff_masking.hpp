//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/read_stream_vector.hpp"
#include "modules/graph_construction.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "modules/graph_read_correction.hpp"
#include "test_utils.hpp"

#include "coloring.hpp"
#include "colored_graph_construction.hpp"
#include "gene_analysis.hpp"
#include "repeat_masking.hpp"
#include "visualization.hpp"

//todo deprecated
namespace cap {

inline void SaveAll(ContigStreams& streams, const vector<string>& suffixes,
        const string& out_root) {
    make_dir(out_root);

    streams.reset();
    for (size_t i = 0; i < streams.size(); ++i) {
        if (!suffixes[i].empty()) {
            string output_filename = out_root + suffixes[i];
            auto rc_wrapped = io::RCWrap<Contig>(streams.ptr_at(i));
            io::osequencestream ostream(output_filename);
            Transfer(*rc_wrapped, ostream);
        }
    }
}

//todo changes the graph!!! color edge splitting!!!
template<class gp_t>
void MakeSaves(gp_t& gp, ContigStreams streams, const string& root,
        const vector<string>& suffixes, bool optional = true) {

    SaveAll(streams, suffixes, root);

    static bool make_optional_saves = true;

    if (!make_optional_saves && optional)
        return;

    make_dir(root);

    streams.reset();

    ColorHandler<Graph> coloring(gp.g, streams.size());
    CoordinatesHandler<Graph> coordinates_handler;
    SplitAndColorGraph(gp, coloring, streams);
    FillPositions(gp, streams, coordinates_handler);

    PrintColoredGraphWithColorFilter(gp.g, coloring, gp.edge_pos,
            root + "colored_split_graph");
}

template<class gp_t>
void RefineGP(gp_t& gp, size_t delta = 5) {
    using namespace debruijn_graph;
    INFO("Constructing graph pack for refinement");

    //outdated
    //debruijn_config::simplification::bulge_remover br_config;
    //br_config.max_bulge_length_coefficient = 3;
    //br_config.max_coverage = 1000.;
    //br_config.max_relative_coverage = 1.2;
    //br_config.max_delta = double(delta);
    //br_config.max_relative_delta = 0.1;

    INFO("Removing bulges");
    debruijn::simplification::RemoveBulges(gp.g, br_config);

    INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");

    debruijn_config::simplification::complex_bulge_remover cbr_config;
    cbr_config.enabled = true;
    cbr_config.pics_enabled = false;
    cbr_config.folder = "";
    cbr_config.max_relative_length = 3;
    cbr_config.max_length_difference = delta;

    INFO("Removing complex bulges");
    debruijn::simplification::RemoveComplexBulges(gp.g, cbr_config);

    TipsProjector<gp_t> tip_projector(gp);
    boost::function<void(EdgeId)> projecting_callback = boost::bind(
            &TipsProjector<gp_t>::ProjectTip, &tip_projector, _1);
    debruijn_config::simplification::tip_clipper tc_config;

    tc_config.condition = "{ tc_lb 2. }";

    INFO("Clipping tips with projection");

    debruijn::simplification::SimplifInfoContainer info_container;

    debruijn::simplification::ClipTipsWithProjection(gp, tc_config, info_container);

    INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");
}

template<class gp_t>
void ConstructGPForRefinement(gp_t& gp, ContigStreams& contigs,
        size_t delta = 5) {
    using namespace debruijn_graph;
    typedef typename gp_t::graph_t Graph;
    INFO("Constructing graph pack for refinement");

    CapConstructGraph(gp.k_value, contigs, gp.g, gp.index);

    RefineGP(gp, delta);
}

template<class gp_t>
ContigStreams RefinedStreams(ContigStreams& streams, const gp_t& gp) {
    ContigStreams refined_streams;
    for (size_t i = 0; i < streams.size(); ++i) {
        refined_streams.push_back(
                make_shared<io::ModifyingWrapper<Contig>>(
                        streams.ptr_at(i),
                        GraphReadCorrectorInstance(gp.g, *MapperInstance(gp))));
    }
    return refined_streams;
}

template<class Seq>
ContigStreams RefineStreams(ContigStreams& streams,
                               size_t k,
                               size_t delta = 5,
                               const std::string &workdir = "tmp") {
    typedef debruijn_graph::KmerStoringEdgeIndex<Graph, Seq, kmer_index_traits<RtSeq>, debruijn_graph::SimpleStoring> RefiningIndex;
    typedef graph_pack<ConjugateDeBruijnGraph, Seq, RefiningIndex> refining_gp_t;
    refining_gp_t gp(k, workdir);

    CapConstructGraph(gp.k_value, streams, gp.g, gp.index);

    RefineGP(gp, delta);

    return RefineStreams(streams, gp);

}


template<class Seq>
void RefineData(const string& base_path,
                            const vector<string>& suffixes,
                            const string& out_root,
                            size_t k,
                            size_t delta = 5,
                            const std::string &workdir = "tmp") {
    ContigStreams streams = OpenStreams(base_path, suffixes, true);
    ContigStreams refined = RefineStreams<Seq>(streams, k, delta, workdir);
    SaveAll(refined, suffixes, out_root);
}

//template<class gp_t>
//void ConstructGPForRefinement(gp_t& gp,
//        io::IReader<io::SingleRead>& raw_stream_1,
//        io::IReader<io::SingleRead>& raw_stream_2, size_t delta = 5) {
//    ContigStreamsPtr streams_ptr = make_shared<ContigStreams>(
//            vector<ContigStream*> { &raw_stream_1, &raw_stream_2 }, false);
//    ConstructGPForRefinement(gp, streams_ptr, delta);
//}

//template<size_t k, class Seq>
//pair<Sequence, Sequence> CorrectGenomes(const Sequence& genome1,
//        const Sequence& genome2, size_t delta = 5) {
//    io::VectorReader<io::SingleRead> stream1(
//            io::SingleRead("first", genome1.str()));
//    io::VectorReader<io::SingleRead> stream2(
//            io::SingleRead("second", genome2.str()));
//
//    typedef debruijn_graph::graph_pack<debruijn_graph::ConjugateDeBruijnGraph,
//            Seq> refining_gp_t;
//    refining_gp_t refining_gp(k, "tmp");
//    ConstructGPForRefinement(refining_gp, stream1, stream2, delta);
//
//    io::ModifyingWrapper<io::SingleRead> refined_stream1(stream1,
//            GraphReadCorrectorInstance(refining_gp.g,
//                    *MapperInstance(refining_gp)));
//    io::ModifyingWrapper<io::SingleRead> refined_stream2(stream2,
//            GraphReadCorrectorInstance(refining_gp.g,
//                    *MapperInstance(refining_gp)));
//
//    pair<Sequence, Sequence> answer = make_pair(FirstSequence(refined_stream1),
//            FirstSequence(refined_stream2));
//    return answer;
//}

//template<size_t k>
//pair<Sequence, Sequence> CorrectGenomes(const pair<Sequence, Sequence>& genomes,
//        size_t delta = 5) {
//    return CorrectGenomes<k>(genomes.first, genomes.second, delta);
//}

//template<size_t k, class Seq>
//pair<Sequence, vector<Sequence>> RefineData(
//        const pair<Sequence, vector<Sequence>>& data) {
//    io::VectorReader<io::SingleRead> stream1(
//            io::SingleRead("first", data.first.str()));
//    io::VectorReader<io::SingleRead> stream2(MakeReads(data.second));
//
//    typedef graph_pack<ConjugateDeBruijnGraph, Seq> refining_gp_t;
//    refining_gp_t refining_gp(k, "tmp");
//    ConstructGPForRefinement(refining_gp, stream1, stream2);
//
//    io::ModifyingWrapper<io::SingleRead> refined_stream1(stream1,
//            GraphReadCorrectorInstance(refining_gp.g,
//                    *MapperInstance(refining_gp)));
//    io::ModifyingWrapper<io::SingleRead> refined_stream2(stream2,
//            GraphReadCorrectorInstance(refining_gp.g,
//                    *MapperInstance(refining_gp)));
//
//    return make_pair(FirstSequence(refined_stream1),
//            AllSequences(refined_stream2));
//}

template<class Seq>
void PerformRefinement(ContigStreams& streams, const string& root,
        const vector<string>& suffixes, size_t k, const string& gene_root,
        GeneCollection& gene_collection) {
    VERIFY(streams.size() == suffixes.size());

    const size_t delta = std::max(size_t(5), k /*/ 5*/);
    typedef graph_pack<ConjugateDeBruijnGraph, Seq, KmerStoringEdgeIndex<Graph, Seq, kmer_index_traits<Seq>, SimpleStoring>> gp_t;
    typedef NewExtendedSequenceMapper<Graph, typename gp_t::index_t> Mapper;

    make_dir(root);
    INFO("Constructing graph pack for k=" << k << " delta=" << delta);
    gp_t gp(unsigned(k), "tmp", 0);

    CapConstructGraph(streams, gp.g, gp.index);

    MakeSaves(gp, streams, root + "before_refinement/", suffixes);

    RefineGP(gp, delta);

    ContigStreams refined_streams = RefinedStreams(streams, gp);

    MakeSaves(gp, refined_streams, root + "after_refinement/", suffixes);

    //todo temporary
    if (!gene_root.empty()) {
        gene_collection.Update(gp);
        ColorHandler<typename gp_t::graph_t> coloring(gp.g);
        string gene_save_dir = root + "updated_gene_info/";
        make_dir(gene_save_dir);
        gene_collection.Save(gene_save_dir, "genomes/", "gene_info.txt");
        string gene_pics_dir = gene_save_dir + "pics/";
        make_dir(gene_pics_dir);
//        WriteGeneLocality(gene_collection, gp, gene_pics_dir, coloring);
    }
    //end temporary
}

inline void PerformIterativeRefinement(ContigStreams& streams,
        const vector<string>& suffixes, const string& out_root,
        vector<size_t> &k_values, const string& gene_root,
        GeneCollection& gene_collection) {

    if (k_values.size() == 0) {
        SaveAll(streams, suffixes, out_root + "final_contigs/");
        return;
    }

    size_t current_k = k_values.back();
    k_values.pop_back();

    string root = out_root + std::to_string(current_k) + "/";

    if (utils::NeedToUseLongSeq(current_k)) {
        omp_set_num_threads(1);
        PerformRefinement<LSeq>(streams, root, suffixes, current_k, gene_root,
                gene_collection);
    } else {
        omp_set_num_threads(8);
        PerformRefinement<RtSeq>(streams, root, suffixes, current_k,
                gene_root, gene_collection);
    }

    ContigStreams corr_streams = OpenStreams(root + "after_refinement/",
            suffixes, true);
    //recursive call
    GeneCollection updated_collection;

    if (!gene_root.empty()) {
        updated_collection.Load(gene_root + "genome_list.txt",
                root + "updated_gene_info/genomes/",
                root + "updated_gene_info/gene_info.txt",
                gene_root + "interesting_orthologs.txt");
    }
    PerformIterativeRefinement(corr_streams, suffixes, out_root, k_values,
            gene_root, updated_collection);

}

inline void PerformIterativeRefinement(const string& base_path,
        const vector<string>& suffixes, const string& out_root,
        vector<size_t>& k_values, bool /* gene_analysis  */= false) {
//    remove_dir(out_root);
    utils::MakeDirPath(out_root);
    ContigStreams streams = OpenStreams(base_path, suffixes, true);

    //stab
    GeneCollection gene_collection;
    PerformIterativeRefinement(streams, suffixes, out_root, k_values, "", gene_collection);
}

//todo temporary
inline void PerformIterativeGeneAnalysis(const string& base_path,
        const string& out_root,
        vector<size_t>& k_values) {

    GeneCollection gene_collection;
    gene_collection.Load(base_path + "genome_list.txt", base_path + "/genomes/",
            base_path + "gene_info.txt",
            base_path + "interesting_orthologs.txt");
    ContigStreams streams;
    vector<string> suffixes;
    for (auto it = gene_collection.genomes.begin(); it != gene_collection.genomes.end(); ++it) {
        streams.push_back(make_shared<io::VectorReadStream<Contig>>(Contig(it->second.name, it->second.sequence.str())));
        suffixes.push_back(it->second.name);
    }
    PerformIterativeRefinement(streams, suffixes, out_root, k_values, base_path,
            gene_collection);
}

}
