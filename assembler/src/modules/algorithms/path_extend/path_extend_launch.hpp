//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * lc_launch.hpp
 *
 *  Created on: Dec 1, 2011
 *      Author: andrey
 */

#ifndef PATH_EXTEND_LAUNCH_HPP_
#define PATH_EXTEND_LAUNCH_HPP_

#include "scaffolder2015/scaffold_graph_constructor.hpp"
#include "pe_config_struct.hpp"
#include "pe_resolver.hpp"
#include "path_extender.hpp"
#include "pe_io.hpp"
#include "path_visualizer.hpp"
#include "loop_traverser.hpp"
#include "assembly_graph/graph_alignment/long_read_storage.hpp"
#include "next_path_searcher.hpp"
#include "scaffolder2015/extension_chooser2015.hpp"
#include "algorithms/genome_consistance_checker.hpp"
#include "scaffolder2015/scaffold_graph.hpp"
#include "scaffolder2015/scaffold_graph_visualizer.hpp"

namespace path_extend {

using namespace debruijn_graph;
typedef omnigraph::de::PairedInfoIndicesT<Graph> PairedInfoIndicesT;

inline size_t FindMaxOverlapedLen(const vector<shared_ptr<PairedInfoLibrary> >& libs) {
    size_t max = 0;
    for (size_t i = 0; i < libs.size(); ++i) {
        max = std::max(libs[i]->GetISMax(), max);
    }
    return max;
}

inline string GetEtcDir(const std::string& output_dir) {
    return output_dir + cfg::get().pe_params.etc_dir + "/";
}

inline void DebugOutputPaths(const conj_graph_pack& gp,
                      const std::string& output_dir, const PathContainer& paths,
                      const string& name) {
    PathInfoWriter path_writer;
    PathVisualizer visualizer;

    DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(gp.g);
    DefaultContigConstructor<ConjugateDeBruijnGraph> constructor(gp.g, corrector);
    ContigWriter writer(gp.g, constructor);

    string etcDir = GetEtcDir(output_dir);
    if (!cfg::get().pe_params.debug_output) {
        return;
    }
    writer.OutputPaths(paths, etcDir + name);
    if (cfg::get().pe_params.output.write_paths) {
        path_writer.WritePaths(paths, etcDir + name + ".dat");
    }
    if (cfg::get().pe_params.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + name + ".dot", name,
                                             paths);
    }
}

inline double GetWeightThreshold(shared_ptr<PairedInfoLibrary> lib, const pe_config::ParamSetT& pset) {
    return lib->IsMp() ? pset.mate_pair_options.weight_threshold : pset.extension_options.weight_threshold;
}

inline double GetPriorityCoeff(shared_ptr<PairedInfoLibrary> lib, const pe_config::ParamSetT& pset) {
    return lib->IsMp() ? pset.mate_pair_options.priority_coeff : pset.extension_options.priority_coeff;
}

inline void SetSingleThresholdForLib(shared_ptr<PairedInfoLibrary> lib, const pe_config::ParamSetT &pset, double threshold, double correction_coeff = 1.0) {
    if  (lib->IsMp()) {
        lib->SetSingleThreshold(pset.mate_pair_options.use_default_single_threshold || math::le(threshold, 0.0) ?
                                pset.mate_pair_options.single_threshold : threshold);
    }
    else {
        double t = pset.extension_options.use_default_single_threshold || math::le(threshold, 0.0) ?
                   pset.extension_options.single_threshold : threshold;
        t = correction_coeff * t;
        lib->SetSingleThreshold(t);
    }
}


inline string MakeNewName(const std::string& contigs_name, const std::string& subname) {
    return contigs_name.substr(0, contigs_name.rfind(".fasta")) + "_" + subname + ".fasta";
}

inline void OutputBrokenScaffolds(PathContainer& paths, int k,
                           const ContigWriter& writer,
                           const std::string& filename) {
    if (!cfg::get().pe_params.param_set.scaffolder_options.on
            or !cfg::get().use_scaffolder
            or cfg::get().pe_params.obs == obs_none) {
        return;
    }

    int min_gap = cfg::get().pe_params.obs == obs_break_all ? k / 2 : k;

    ScaffoldBreaker breaker(min_gap);
    breaker.Split(paths);
    breaker.container().SortByLength();
    writer.OutputPaths(breaker.container(), filename);
}

inline void AddPathsToContainer(const conj_graph_pack& gp,
                         const std::vector<PathInfo<Graph> > paths,
                         size_t size_threshold, PathContainer& result) {
    for (size_t i = 0; i < paths.size(); ++i) {
        auto path = paths.at(i);
        vector<EdgeId> edges = path.getPath();
        if (edges.size() <= size_threshold) {
            continue;
        }
        BidirectionalPath* new_path = new BidirectionalPath(gp.g, edges);
        BidirectionalPath* conj_path = new BidirectionalPath(new_path->Conjugate());
        new_path->SetWeight((float) path.getWeight());
        conj_path->SetWeight((float) path.getWeight());
        result.AddPair(new_path, conj_path);
    }
    DEBUG("Long reads paths " << result.size() << " == ");
}

double GetSingleReadsFilteringThreshold(const io::LibraryType& type) {
    if (type == io::LibraryType::PacBioReads || type == io::LibraryType::SangerReads || type == io::LibraryType::NanoporeReads) {
        return cfg::get().pe_params.long_reads.pacbio_reads.filtering;
    }  else if (type == io::LibraryType::PathExtendContigs){
        return cfg::get().pe_params.long_reads.meta_contigs.filtering;
    } else if (io::SequencingLibraryBase::is_contig_lib(type)) {
        return cfg::get().pe_params.long_reads.contigs.filtering;
    }
    return cfg::get().pe_params.long_reads.single_reads.filtering;
}

double GetSingleReadsWeightPriorityThreshold(const io::LibraryType& type) {
    if (type == io::LibraryType::PacBioReads || type == io::LibraryType::SangerReads || type == io::LibraryType::NanoporeReads) {
        return cfg::get().pe_params.long_reads.pacbio_reads.weight_priority;
    } else if (type == io::LibraryType::PathExtendContigs){
        return cfg::get().pe_params.long_reads.meta_contigs.weight_priority;
    } else if (io::SequencingLibraryBase::is_contig_lib(type)) {
        return cfg::get().pe_params.long_reads.contigs.weight_priority;
    }
    return cfg::get().pe_params.long_reads.single_reads.weight_priority;
}

double GetSingleReadsUniqueEdgePriorityThreshold(const io::LibraryType& type) {
    if (cfg::get().ds.single_cell &&
            (type == io::LibraryType::PacBioReads || type == io::LibraryType::SangerReads || type == io::LibraryType::NanoporeReads)) {
        return 10000.0;
    }
    if (type == io::LibraryType::PacBioReads || type == io::LibraryType::SangerReads || type == io::LibraryType::NanoporeReads) {
        return cfg::get().pe_params.long_reads.pacbio_reads.unique_edge_priority;
    } else if (io::SequencingLibraryBase::is_contig_lib(type)) {
        return cfg::get().pe_params.long_reads.contigs.unique_edge_priority;
    }
    return cfg::get().pe_params.long_reads.single_reads.unique_edge_priority;
}

bool HasOnlyMPLibs() {
    for (const auto& lib : cfg::get().ds.reads) {
        if (!((lib.type() == io::LibraryType::MatePairs || lib.type() == io::LibraryType::HQMatePairs) &&
              lib.data().mean_insert_size > 0.0)) {
            return false;
        }
    }
    return true;
}

bool UseCoverageResolverForSingleReads(const io::LibraryType& type) {
    return HasOnlyMPLibs() && (type == io::LibraryType::HQMatePairs);
}

inline size_t CountEdgesInGraph(const Graph& g) {
    size_t count = 0;
    for (auto iter = g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        count++;
    }
    return count;
}

inline size_t GetNumberMPPaths(const Graph& g) {
    size_t count_edge = CountEdgesInGraph(g);
    if (count_edge < 1000) {
        return 1000;
    }
    if (count_edge < 10000) {
        return 100;
    }
    return 50;
}

inline string LibStr(size_t count) {
    return count == 1 ? "library" : "libraries";
}

inline void ClonePathContainer(PathContainer& spaths, PathContainer& tpaths, GraphCoverageMap& tmap) {
    tpaths.clear();
    tmap.Clear();

    for (auto iter = spaths.begin(); iter != spaths.end(); ++iter) {
        BidirectionalPath& path = *iter.get();
        BidirectionalPath* new_path = new BidirectionalPath(path.graph());
        new_path->Subscribe(&tmap);
        new_path->PushBack(path);

        BidirectionalPath& cpath = *iter.getConjugate();
        BidirectionalPath* new_cpath = new BidirectionalPath(cpath.graph());
        new_cpath->Subscribe(&tmap);
        new_cpath->PushBack(cpath);

        tpaths.AddPair(new_path, new_cpath);
    }
}

inline void FinalizePaths(PathContainer& paths, GraphCoverageMap& cover_map, size_t min_edge_len, size_t max_path_diff, bool mate_pairs = false) {
    PathExtendResolver resolver(cover_map.graph());


    if (cfg::get().pe_params.param_set.remove_overlaps) {
        resolver.removeOverlaps(paths, cover_map, min_edge_len, max_path_diff, cfg::get().pe_params.param_set.cut_all_overlaps);
    }
    else {
        resolver.removeEqualPaths(paths, cover_map, min_edge_len);
    }
    if (mate_pairs) {
        resolver.RemoveMatePairEnds(paths, min_edge_len);
    }
    if (cfg::get().avoid_rc_connections) {
        paths.FilterInterstandBulges();
    }
    paths.FilterEmptyPaths();
    if (!mate_pairs) {
        resolver.addUncoveredEdges(paths, cover_map);
    }
    paths.SortByLength();
    for(auto& path : paths) {
        path.first->ResetOverlaps();
    }

}

inline void TraverseLoops(PathContainer& paths, GraphCoverageMap& cover_map, shared_ptr<ContigsMaker> extender) {
    INFO("Traversing tandem repeats");
    LoopTraverser loopTraverser(cover_map.graph(), cover_map, extender);
    loopTraverser.TraverseAllLoops();
    paths.SortByLength();
}

inline bool IsForSingleReadExtender(const io::SequencingLibrary<debruijn_config::DataSetData> &lib) {
    io::LibraryType lt = lib.type();
    return (lib.data().single_reads_mapped ||
            lt == io::LibraryType::PacBioReads ||
            lt == io::LibraryType::SangerReads ||
            lt == io::LibraryType::NanoporeReads ||
            lib.is_contig_lib());
}

inline bool IsForPEExtender(const io::SequencingLibrary<debruijn_config::DataSetData> &lib) {
    return (lib.type() == io::LibraryType::PairedEnd &&
            lib.data().mean_insert_size > 0.0);
}

inline bool IsForShortLoopExtender(const io::SequencingLibrary<debruijn_config::DataSetData> &lib) {
    return (lib.type() == io::LibraryType::PairedEnd &&
            lib.data().mean_insert_size > 0.0);
}

inline bool IsForScaffoldingExtender(const io::SequencingLibrary<debruijn_config::DataSetData> &lib) {
    return (lib.type() == io::LibraryType::PairedEnd &&
            lib.data().mean_insert_size > 0.0);
}

inline bool IsForMPExtender(const io::SequencingLibrary<debruijn_config::DataSetData> &lib) {
    return lib.data().mean_insert_size > 0.0 &&
            (lib.type() == io::LibraryType::HQMatePairs ||
             lib.type() == io::LibraryType::MatePairs);
}

enum class PathExtendStage {
    PEStage,
    PEPolishing,
    MPStage,
    FinalizingPEStage,
    FinalPolishing,
};

inline bool IsPEStage(PathExtendStage stage) {
    return stage == PathExtendStage::PEPolishing || stage == PathExtendStage::PEStage;
}

inline bool IsMPStage(PathExtendStage stage) {
    return stage == PathExtendStage::MPStage;
}

inline bool IsFinalStage(PathExtendStage stage) {
    return stage == PathExtendStage::FinalizingPEStage || stage == PathExtendStage::FinalPolishing;
}

inline bool IsPolishingStage(PathExtendStage stage) {
    return stage == PathExtendStage::PEPolishing || stage == PathExtendStage::FinalPolishing;
}


template<class Index>
inline shared_ptr<PairedInfoLibrary> MakeNewLib(const conj_graph_pack::graph_t& g,
                                     const Index& paired_index,
                                     size_t index) {
    const auto& lib = cfg::get().ds.reads[index];
    size_t read_length = lib.data().read_length;
    size_t is = (size_t) lib.data().mean_insert_size;
    int is_min = (int) lib.data().insert_size_left_quantile;
    int is_max = (int) lib.data().insert_size_right_quantile;
    int var = (int) lib.data().insert_size_deviation;
    bool is_mp = lib.type() == io::LibraryType::MatePairs ||  lib.type() == io::LibraryType::HQMatePairs ;
    return make_shared< PairedInfoLibraryWithIndex<decltype(paired_index[index])> >(cfg::get().K, g, read_length,
                                                                                    is, is_min > 0.0 ? size_t(is_min) : 0, is_max > 0.0 ? size_t(is_max) : 0,
                                                                                    size_t(var),
                                                                                    paired_index[index], is_mp,
                                                                                    lib.data().insert_size_distribution);
}

inline shared_ptr<SimpleExtender> MakeLongReadsExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, size_t lib_index,
                                                        const pe_config::ParamSetT& pset) {
    PathContainer paths;
    AddPathsToContainer(gp, gp.single_long_reads[lib_index].GetAllPaths(), 1, paths);

    const auto& lib = cfg::get().ds.reads[lib_index];
    shared_ptr<ExtensionChooser> longReadEC =
            make_shared<LongReadsExtensionChooser>(gp.g, paths, GetSingleReadsFilteringThreshold(lib.type()),
                                                   GetSingleReadsWeightPriorityThreshold(lib.type()),
                                                   GetSingleReadsUniqueEdgePriorityThreshold(lib.type()),
                                                   pset.extension_options.max_repeat_length);

    size_t resolvable_repeat_length_bound = 10000ul;
    if (!lib.is_contig_lib()) {
        resolvable_repeat_length_bound = std::max(resolvable_repeat_length_bound, lib.data().read_length);
    }
    INFO("resolvable_repeat_length_bound set to " << resolvable_repeat_length_bound);
    return make_shared<SimpleExtender>(gp, cov_map, longReadEC, resolvable_repeat_length_bound,  
            pset.loop_removal.max_loops, true, UseCoverageResolverForSingleReads(lib.type()));
}

inline shared_ptr<SimpleExtender> MakeLongEdgePEExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                                         size_t lib_index, const pe_config::ParamSetT& pset, bool investigate_loops) {
    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.clustered_indices, lib_index);
    SetSingleThresholdForLib(lib, pset, cfg::get().ds.reads[lib_index].data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << lib->GetSingleThreshold());

    shared_ptr<WeightCounter> wc = make_shared<PathCoverWeightCounter>(gp.g, lib, pset.normalize_weight);
    shared_ptr<ExtensionChooser> extension = make_shared<LongEdgeExtensionChooser>(gp.g, wc, GetWeightThreshold(lib, pset), GetPriorityCoeff(lib, pset));
    return make_shared<SimpleExtender>(gp, cov_map, extension, lib->GetISMax(), pset.loop_removal.max_loops, investigate_loops, false);
}


inline shared_ptr<SimpleExtender> MakeMetaExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                       size_t lib_index, const pe_config::ParamSetT& pset, bool investigate_loops) {
    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.clustered_indices, lib_index);
    VERIFY(!lib->IsMp());

    shared_ptr<WeightCounter> wc = make_shared<MetagenomicWeightCounter>(gp.g, lib, /*read_length*/cfg::get().ds.RL(), 
                            /*normalized_threshold*/ 0.3, /*raw_threshold*/ 3, /*estimation_edge_length*/ 300);
    shared_ptr<SimpleExtensionChooser> extension = make_shared<SimpleExtensionChooser>(gp.g, wc, 
                                                        pset.extension_options.weight_threshold, 
                                                        pset.extension_options.priority_coeff);
    return make_shared<SimpleExtender>(gp, cov_map, extension, lib->GetISMax(), pset.loop_removal.max_loops, investigate_loops, false);
}

inline shared_ptr<SimpleExtender> MakePEExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                       size_t lib_index, const pe_config::ParamSetT& pset, bool investigate_loops) {
    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.clustered_indices, lib_index);
    SetSingleThresholdForLib(lib, pset, cfg::get().ds.reads[lib_index].data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << lib->GetSingleThreshold());

    shared_ptr<WeightCounter> wc = make_shared<PathCoverWeightCounter>(gp.g, lib, pset.normalize_weight);
    shared_ptr<SimpleExtensionChooser> extension = make_shared<SimpleExtensionChooser>(gp.g, wc, GetWeightThreshold(lib, pset), GetPriorityCoeff(lib, pset));
    return make_shared<SimpleExtender>(gp, cov_map, extension, lib->GetISMax(), pset.loop_removal.max_loops, investigate_loops, false);
}

inline shared_ptr<PathExtender> MakeScaffoldingExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                       size_t lib_index, const pe_config::ParamSetT& pset) {
    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.scaffolding_indices, lib_index);

    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(gp.g, lib);
    //FIXME this variable was not used!
    //double prior_coef = GetPriorityCoeff(lib, pset);
    //FIXME review parameters
    //todo put parameters in config
    //FIXME remove max_must_overlap from config
    double var_coeff = 3.0;
    auto scaff_chooser = std::make_shared<ScaffoldingExtensionChooser>(gp.g, counter, var_coeff);

    vector<shared_ptr<GapJoiner>> joiners;

    if (pset.scaffolder_options.use_la_gap_joiner) {
        joiners.push_back(std::make_shared<LAGapJoiner>(gp.g, pset.scaffolder_options.min_overlap_length,
                                                    pset.scaffolder_options.flank_multiplication_coefficient,
                                                    pset.scaffolder_options.flank_addition_coefficient));
    }

    joiners.push_back(std::make_shared<HammingGapJoiner>(gp.g, pset.scaffolder_options.min_gap_score,
                                                 pset.scaffolder_options.short_overlap,
                                                 (int) 2 * cfg::get().ds.RL()));

    auto composite_gap_joiner = std::make_shared<CompositeGapJoiner>(gp.g, 
                                                joiners, 
                                                size_t(pset.scaffolder_options.max_can_overlap * (double) gp.g.k()),
                                                int(math::round((double) gp.g.k() - var_coeff * (double) lib->GetIsVar())),
                                                pset.scaffolder_options.artificial_gap);

    return make_shared<ScaffoldingPathExtender>(gp, cov_map, scaff_chooser, composite_gap_joiner, lib->GetISMax(), pset.loop_removal.max_loops, false);
}


inline shared_ptr<PathExtender> MakeScaffolding2015Extender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                                        size_t lib_index, const pe_config::ParamSetT& pset, const ScaffoldingUniqueEdgeStorage& storage) {
    shared_ptr<PairedInfoLibrary> lib;
    INFO("for lib " << lib_index);

    //TODO:: temporary solution
    if (gp.paired_indices[lib_index].size() > gp.clustered_indices[lib_index].size()) {
        INFO("Paired unclustered indices not empty, using them");
        lib = MakeNewLib(gp.g, gp.paired_indices, lib_index);
    } else if (gp.clustered_indices[lib_index].size() != 0 ) {
        INFO("clustered indices not empty, using them");
        lib = MakeNewLib(gp.g, gp.clustered_indices, lib_index);
    } else {
        ERROR("All paired indices are empty!");
    }

    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(gp.g, lib);
//TODO::was copypasted from MakeScaffoldingExtender
//TODO::REWRITE
    double var_coeff = 3.0;
    DEBUG("here creating extchooser");
//TODO: 2 is relative weight cutoff, to config!
    auto scaff_chooser = std::make_shared<ExtensionChooser2015>(gp.g, counter, var_coeff, storage, 2, lib_index);

    auto gap_joiner = std::make_shared<HammingGapJoiner>(gp.g, pset.scaffolder_options.min_gap_score,
                                                         pset.scaffolder_options.short_overlap,
                                                         (int) 2 * cfg::get().ds.RL());

    return make_shared<ScaffoldingPathExtender>(gp, cov_map, scaff_chooser, gap_joiner, lib->GetISMax(), pset.loop_removal.max_loops, false , false);
}


inline shared_ptr<SimpleExtender> MakeMPExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, const PathContainer& paths,
                                       size_t lib_index, const pe_config::ParamSetT& pset) {

    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.paired_indices, lib_index);
    SetSingleThresholdForLib(lib, pset, cfg::get().ds.reads[lib_index].data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << lib->GetSingleThreshold());

    size_t max_number_of_paths_to_search = GetNumberMPPaths(gp.g);
    DEBUG("max number of mp paths " << max_number_of_paths_to_search);

    shared_ptr<MatePairExtensionChooser> chooser = make_shared<MatePairExtensionChooser>(gp.g, lib, paths, max_number_of_paths_to_search);
    return make_shared<SimpleExtender>(gp, cov_map, chooser, lib->GetISMax(), pset.loop_removal.mp_max_loops, true, false);
}

inline shared_ptr<SimpleExtender> MakeCoordCoverageExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                       const pe_config::ParamSetT& pset) {
    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.paired_indices, 0);
    CoverageAwareIdealInfoProvider provider(gp.g, lib, 1000, 2000);
    shared_ptr<CoordinatedCoverageExtensionChooser> chooser = make_shared<CoordinatedCoverageExtensionChooser>(gp.g, provider,
            pset.coordinated_coverage.max_edge_length_in_repeat, pset.coordinated_coverage.delta);
    return make_shared<SimpleExtender>(gp, cov_map, chooser, -1ul, pset.loop_removal.mp_max_loops, true, false);
}

inline shared_ptr<SimpleExtender> MakeRNAExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                                 size_t lib_index, const pe_config::ParamSetT& pset, bool investigate_loops) {
    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.clustered_indices, lib_index);
    SetSingleThresholdForLib(lib, pset, cfg::get().ds.reads[lib_index].data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << lib->GetSingleThreshold());

    shared_ptr<WeightCounter> wc = make_shared<PathCoverWeightCounter>(gp.g, lib, pset.normalize_weight);
    shared_ptr<RNAExtensionChooser> extension = make_shared<RNAExtensionChooser>(gp.g, wc, GetWeightThreshold(lib, pset), GetPriorityCoeff(lib, pset));
    return make_shared<MultiExtender>(gp, cov_map, extension, lib->GetISMax(), pset.loop_removal.max_loops, investigate_loops, false);
}

inline shared_ptr<SimpleExtender> MakeRNALongReadsExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, size_t lib_index,
                                                        const pe_config::ParamSetT& pset) {
    VERIFY_MSG(false, "Long reads rna extender us not implemented yet")
    PathContainer paths;
    AddPathsToContainer(gp, gp.single_long_reads[lib_index].GetAllPaths(), 1, paths);

    const auto& lib = cfg::get().ds.reads[lib_index];
    shared_ptr<ExtensionChooser> longReadEC =
        make_shared<LongReadsExtensionChooser>(gp.g, paths, GetSingleReadsFilteringThreshold(lib.type()),
                                               GetSingleReadsWeightPriorityThreshold(lib.type()),
                                               GetSingleReadsUniqueEdgePriorityThreshold(lib.type()),
                                               pset.extension_options.max_repeat_length);

    size_t resolvable_repeat_length_bound = 10000ul;
    if (!lib.is_contig_lib()) {
        resolvable_repeat_length_bound = std::max(resolvable_repeat_length_bound, lib.data().read_length);
    }
    INFO("resolvable_repeat_length_bound set to " << resolvable_repeat_length_bound);
    return make_shared<SimpleExtender>(gp, cov_map, longReadEC, resolvable_repeat_length_bound,
                                       pset.loop_removal.max_loops, true, UseCoverageResolverForSingleReads(lib.type()));
}



inline bool InsertSizeCompare(const shared_ptr<PairedInfoLibrary> lib1,
                              const shared_ptr<PairedInfoLibrary> lib2) {
    return lib1->GetISMax() < lib2->GetISMax();
}

template<typename Base, typename T>
inline bool instanceof(const T *ptr) {
    return dynamic_cast<const Base*>(ptr) != nullptr;
}

//Used for debug purpose only
inline void PrintExtenders(vector<shared_ptr<PathExtender> >& extenders) {
    DEBUG("Extenders in vector:");
    for(size_t i = 0; i < extenders.size(); ++i) {
        string type = typeid(*extenders[i]).name();
        DEBUG("Extender #i" << type);
        if (instanceof<SimpleExtender>(extenders[i].get())) {
            auto ec = ((SimpleExtender *) extenders[i].get())->GetExtensionChooser();
            string chooser_type = typeid(*ec).name();
            DEBUG("    Extender #i" << chooser_type);
        }
        else if (instanceof<ScaffoldingPathExtender>(extenders[i].get())) {
            auto ec = ((ScaffoldingPathExtender *) extenders[i].get())->GetExtensionChooser();
            string chooser_type = typeid(*ec).name();
            DEBUG("    Extender #i" << chooser_type);
        }
    }
}

inline vector<shared_ptr<PathExtender> > MakeAllExtenders(PathExtendStage stage, const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                            const pe_config::ParamSetT& pset, const ScaffoldingUniqueEdgeStorage& storage, const PathContainer& paths_for_mp = PathContainer()) {

    vector<shared_ptr<PathExtender> > result;
    vector<shared_ptr<PathExtender> > pes;
    vector<shared_ptr<PathExtender> > pes2015;
    vector<shared_ptr<PathExtender> > pe_loops;
    vector<shared_ptr<PathExtender> > pe_scafs;
    vector<shared_ptr<PathExtender> > mps;

    size_t single_read_libs = 0;
    size_t pe_libs = 0;
    size_t scf_pe_libs = 0;
    size_t mp_libs = 0;

    for (io::LibraryType lt : io::LibraryPriotity) {
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            const auto& lib = cfg::get().ds.reads[i];            
            if (lib.type() != lt)
                continue;

            //TODO: scaff2015 does not need any single read libs?
            if (IsForSingleReadExtender(lib) && pset.sm != sm_2015) {
                result.push_back(MakeLongReadsExtender(gp, cov_map, i, pset));
                ++single_read_libs;
            }
            if (IsForPEExtender(lib)) {
                ++pe_libs;
                if (IsPEStage(stage) && (pset.sm == sm_old_pe_2015 || pset.sm == sm_old || pset.sm == sm_combined)) {
                    if (cfg::get().ds.meta)
                        //TODO proper configuration via config
                        pes.push_back(MakeMetaExtender(gp, cov_map, i, pset, false));
                    else if (cfg::get().ds.moleculo)
                        pes.push_back(MakeLongEdgePEExtender(gp, cov_map, i, pset, false));
                    else if (cfg::get().ds.rna && !IsPolishingStage(stage))
                        pes.push_back(MakeRNAExtender(gp, cov_map, i, pset, false));
                    else
                        pes.push_back(MakePEExtender(gp, cov_map, i, pset, false));
                }
                else if (pset.sm == sm_2015) {
                    pes2015.push_back(MakeScaffolding2015Extender(gp, cov_map, i, pset, storage));
                }
            }
            if (IsForShortLoopExtender(lib) && (pset.sm == sm_old_pe_2015 || pset.sm == sm_old || pset.sm == sm_combined)) {
                if (cfg::get().ds.meta)
                    pes.push_back(MakeMetaExtender(gp, cov_map, i, pset, true));
                else if (cfg::get().ds.rna && !IsPolishingStage(stage))
                    pes.push_back(MakeRNAExtender(gp, cov_map, i, pset, true));
                else
                    pe_loops.push_back(MakePEExtender(gp, cov_map, i, pset, true));
            }
            if (IsForScaffoldingExtender(lib) && cfg::get().use_scaffolder && pset.scaffolder_options.on) {
                ++scf_pe_libs;
                if (pset.sm == sm_old || pset.sm == sm_combined) {
                    pe_scafs.push_back(MakeScaffoldingExtender(gp, cov_map, i, pset));
                }
                if (pset.sm == sm_old_pe_2015 || pset.sm == sm_combined) {
                    pe_scafs.push_back(MakeScaffolding2015Extender(gp, cov_map, i, pset, storage));
                }
            }
            if (IsForMPExtender(lib) && IsMPStage(stage)) {
                ++mp_libs;
                if (pset.sm == sm_old || pset.sm == sm_combined) {
                    mps.push_back(MakeMPExtender(gp, cov_map, paths_for_mp, i, pset));
                }
                if (is_2015_scaffolder_enabled(pset.sm)) {
                    mps.push_back(MakeScaffolding2015Extender(gp, cov_map, i, pset, storage));
                }
            }
        }

        //std::sort(scaff_libs.begin(), scaff_libs.end(), InsertSizeCompare);
        result.insert(result.end(), pes.begin(), pes.end());
        result.insert(result.end(), pes2015.begin(), pes2015.end());
        result.insert(result.end(), pe_loops.begin(), pe_loops.end());
        result.insert(result.end(), pe_scafs.begin(), pe_scafs.end());
        result.insert(result.end(), mps.begin(), mps.end());
        pes.clear();
        pe_loops.clear();
        pe_scafs.clear();
        pes2015.clear();
        mps.clear();
    }

    INFO("Using " << pe_libs << " paired-end " << LibStr(pe_libs));
    INFO("Using " << scf_pe_libs << " paired-end scaffolding " << LibStr(scf_pe_libs));
    INFO("Using " << mp_libs << " mate-pair " << LibStr(mp_libs));
    INFO("Using " << single_read_libs << " single read " << LibStr(single_read_libs));
    INFO("Scaffolder is " << (pset.scaffolder_options.on ? "on" : "off"));

    if(pset.use_coordinated_coverage) {
        INFO("Using additional coordinated coverage extender");
        result.push_back(MakeCoordCoverageExtender(gp, cov_map, pset));
    }

    PrintExtenders(result);
    return result;
}

inline shared_ptr<scaffold_graph::ScaffoldGraph> ConstructScaffoldGraph(const conj_graph_pack& gp,
                                                                        const ScaffoldingUniqueEdgeStorage& edge_storage,
                                                                        const pe_config::ParamSetT::ScaffoldGraphParamsT& params) {
    using namespace scaffold_graph;
    vector<shared_ptr<ConnectionCondition>> conditions;

    INFO("Constructing connections");
    LengthEdgeCondition edge_condition(gp.g, edge_storage.GetMinLength());

    for (size_t lib_index = 0; lib_index < cfg::get().ds.reads.lib_count(); ++lib_index) {
        auto lib = cfg::get().ds.reads[lib_index];
        if (lib.is_paired()) {
            shared_ptr<PairedInfoLibrary> paired_lib;
            if (IsForMPExtender(lib))
                paired_lib = MakeNewLib(gp.g, gp.paired_indices, lib_index);
            else if (IsForPEExtender(lib))
                paired_lib = MakeNewLib(gp.g, gp.clustered_indices, lib_index);
            else
                INFO("Unusable paired lib #" << lib_index);
            conditions.push_back(make_shared<AdvancedPairedConnectionCondition>(gp.g, paired_lib, lib_index,
                                                                                params.always_add,
                                                                                params.never_add,
                                                                                params.relative_threshold));
        }
    }
    if (params.graph_connectivity) {
        auto as_con = make_shared<AssemblyGraphConnectionCondition>(gp.g, params.max_path_length, edge_storage);
        for (auto e_iter = gp.g.ConstEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
            if (edge_condition.IsSuitable(*e_iter))
                as_con->AddInterestingEdge(*e_iter);
        }
        conditions.push_back(as_con);
    }
    INFO("Total conditions " << conditions.size());

    INFO("Constructing scaffold graph from set of size " << edge_storage.GetSet().size());

    DefaultScaffoldGraphConstructor constructor(gp.g, edge_storage.GetSet(), conditions, edge_condition);
    auto scaffoldGraph = constructor.Construct();

    INFO("Scaffold graph contains " << scaffoldGraph->VertexCount() << " vertices and " << scaffoldGraph->EdgeCount() << " edges");
    return scaffoldGraph;
}


inline void PrintScaffoldGraph(shared_ptr<scaffold_graph::ScaffoldGraph> scaffoldGraph,
                               const set<EdgeId> main_edge_set,
                               const string& filename) {
    using namespace scaffold_graph;

    auto vcolorer = make_shared<ScaffoldVertexSetColorer>(main_edge_set);
    auto ecolorer = make_shared<ScaffoldEdgeColorer>();
    CompositeGraphColorer <ScaffoldGraph> colorer(vcolorer, ecolorer);

    INFO("Visualizing single grpah");
    ScaffoldGraphVisualizer singleVisualizer(*scaffoldGraph, false);
    std::ofstream single_dot;
    single_dot.open((filename + "_single.dot").c_str());
    singleVisualizer.Visualize(single_dot, colorer);
    single_dot.close();

    INFO("Visualizing paired grpah");
    ScaffoldGraphVisualizer pairedVisualizer(*scaffoldGraph, true);
    std::ofstream paired_dot;
    paired_dot.open((filename + "_paired.dot").c_str());
    pairedVisualizer.Visualize(paired_dot, colorer);
    paired_dot.close();

    INFO("Printing scaffold grpah");
    std::ofstream data_stream;
    data_stream.open((filename + ".data").c_str());
    scaffoldGraph->Print(data_stream);
    data_stream.close();
}


inline size_t FindOverlapLenForStage(PathExtendStage stage) {
    size_t res = 0;
    for (const auto& lib : cfg::get().ds.reads) {
        if (IsForPEExtender(lib) && IsPEStage(stage)) {
            res = max(res, (size_t) lib.data().insert_size_right_quantile);
        } else if (IsForShortLoopExtender(lib)) {
            res = max(res, (size_t) lib.data().insert_size_right_quantile);
        } else if (IsForMPExtender(lib) && IsMPStage(stage)) {
            res = max(res, (size_t) lib.data().insert_size_right_quantile);
        }
    }
    return res;
}

inline bool MPLibsExist() {
    for (const auto& lib : cfg::get().ds.reads)
        if (IsForMPExtender(lib))
            return true;

    return false;
}

inline void CountMisassembliesWithReference(debruijn_graph::GenomeConsistenceChecker& genome_checker, const PathContainer& paths) {
    size_t total_mis = 0 , gap_mis = 0;
    genome_checker.SpellGenome();
    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        BidirectionalPath *path = iter.get();
        auto map_res = genome_checker.CountMisassemblies(*path);
        if (map_res.misassemblies > 0) {
            INFO ("there are " << map_res.misassemblies << " misassemblies in path: ");
            path->PrintInfo();
            total_mis += map_res.misassemblies;
        }
        if (map_res.wrong_gap_size > 0) {
            INFO ("there are " << map_res.wrong_gap_size << " wrong gaps in path: ");
            path->PrintInfo();
            gap_mis += map_res.wrong_gap_size;
        }
    }
    INFO ("In total found " << total_mis << " misassemblies " << " and " << gap_mis << " gaps.");
}

inline ScaffoldingUniqueEdgeStorage FillUniqueEdgeStorage(const conj_graph_pack& gp,
                                                                      size_t& min_unique_length,
                                                                      double& unique_variation) {

    ScaffoldingUniqueEdgeStorage main_unique_storage;
    //Setting scaffolding2015 parameters
    if (cfg::get().pe_params.param_set.scaffolding2015.autodetect) {
        INFO("Autodetecting unique edge set parameters...");
        bool pe_found = false;
        //TODO constants
        size_t min_MP_IS = 10000;
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {

            if (IsForPEExtender(cfg::get().ds.reads[i])) {
                pe_found = true;
            }
            if (IsForMPExtender(cfg::get().ds.reads[i])) {
                min_MP_IS = min(min_MP_IS, (size_t) cfg::get().ds.reads[i].data().mean_insert_size);
            }
        }
        if (pe_found) {
            //TODO constants
            unique_variation = 0.5;
            INFO("PE lib found, we believe in coverage");
        } else {
            unique_variation = 50;
            INFO("No paired libs found, we do not believe in coverage");
        }
        min_unique_length = min_MP_IS;
        INFO("Minimal unique edge length set to the smallest MP library IS: " << min_unique_length);

    } else {
        INFO("Unique edge set constructed with parameters from config : length " << min_unique_length
                 << " variation " << unique_variation);
    }
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp, min_unique_length, unique_variation);
    unique_edge_analyzer.FillUniqueEdgeStorage(main_unique_storage);

    return main_unique_storage;
}

inline void ResolveRepeatsPe(conj_graph_pack& gp,
        const std::string& output_dir,
        const std::string& contigs_name,
        bool traversLoops,
        boost::optional<std::string> broken_contigs) {

    INFO("ExSPAnder repeat resolving tool started");

    ScaffoldingUniqueEdgeStorage main_unique_storage;
    auto sc_mode = cfg::get().pe_params.param_set.sm;
    auto min_unique_length = cfg::get().pe_params.param_set.scaffolding2015.min_unique_length;
    auto unique_variaton = cfg::get().pe_params.param_set.scaffolding2015.unique_coverage_variation;

    if (is_2015_scaffolder_enabled(sc_mode)) {
        main_unique_storage = FillUniqueEdgeStorage(gp, min_unique_length, unique_variaton);
    }

    make_dir(output_dir);
    make_dir(GetEtcDir(output_dir));
    const pe_config::ParamSetT &pset = cfg::get().pe_params.param_set;

    //Scaffold graph
    shared_ptr<scaffold_graph::ScaffoldGraph> scaffoldGraph;
    if (cfg::get().pe_params.param_set.scaffold_graph_params.construct) {
        scaffoldGraph = ConstructScaffoldGraph(gp, main_unique_storage, cfg::get().pe_params.param_set.scaffold_graph_params);
        if (cfg::get().pe_params.param_set.scaffold_graph_params.output) {
            PrintScaffoldGraph(scaffoldGraph, main_unique_storage.GetSet(), GetEtcDir(output_dir) + "scaffold_graph");
        }
    }


    DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(gp.g);
    DefaultContigConstructor<ConjugateDeBruijnGraph> constructor(gp.g, corrector);
    ContigWriter writer(gp.g, constructor, gp.components);


//make pe + long reads extenders
    GraphCoverageMap cover_map(gp.g);
    INFO("SUBSTAGE = paired-end libraries")
    PathExtendStage exspander_stage = PathExtendStage::PEStage;
    vector<shared_ptr<PathExtender> > all_libs = MakeAllExtenders(exspander_stage, gp, cover_map, pset,
                                                                  main_unique_storage);

    //Parameters are subject to change
    size_t max_is_right_quantile = max(FindOverlapLenForStage(exspander_stage), gp.g.k() + 100);
    size_t min_edge_len = 100;

    shared_ptr<CompositeExtender> mainPE = make_shared<CompositeExtender>(gp.g, cover_map, all_libs,
                                                                          max_is_right_quantile, main_unique_storage,
                                                                          cfg::get().pe_params.param_set.extension_options.max_repeat_length);

//extend pe + long reads
    PathExtendResolver resolver(gp.g);
    auto seeds = resolver.makeSimpleSeeds();
    DebugOutputPaths(gp, output_dir, seeds, "init_paths");
    seeds.SortByLength();
    INFO("Growing paths using paired-end and long single reads");
    auto paths = resolver.extendSeeds(seeds, *mainPE);
    paths.SortByLength();
    DebugOutputPaths(gp, output_dir, paths, "pe_before_overlap");

    PathContainer clone_paths;
    GraphCoverageMap clone_map(gp.g);
    bool mp_exist = MPLibsExist();

    if (mp_exist) {
        ClonePathContainer(paths, clone_paths, clone_map);
    }

    exspander_stage = PathExtendStage::PEPolishing;
    all_libs = MakeAllExtenders(exspander_stage, gp, cover_map, pset, main_unique_storage);
    mainPE = make_shared<CompositeExtender>(gp.g, cover_map, all_libs,
                                            max_is_right_quantile, main_unique_storage,
                                            cfg::get().pe_params.param_set.extension_options.max_repeat_length);

    //We do not run overlap removal in 2015 mode
    if (!is_2015_scaffolder_enabled(sc_mode))
        FinalizePaths(paths, cover_map, min_edge_len, max_is_right_quantile);
    if (broken_contigs.is_initialized()) {
        OutputBrokenScaffolds(paths, (int) gp.g.k(), writer,
                              output_dir + (mp_exist ? "pe_contigs" : broken_contigs.get()));
    }
    DebugOutputPaths(gp, output_dir, paths, "pe_before_traverse");
    if (traversLoops) {
        TraverseLoops(paths, cover_map, mainPE);
        FinalizePaths(paths, cover_map, min_edge_len, max_is_right_quantile);
    }
    DebugOutputPaths(gp, output_dir, paths, (mp_exist ? "pe_final_paths" : "final_paths"));
    writer.OutputPaths(paths, output_dir + (mp_exist ? "pe_scaffolds" : contigs_name));

    cover_map.Clear();
    paths.DeleteAllPaths();
    if (!mp_exist) {
        return;
    }

//MP
    DebugOutputPaths(gp, output_dir, clone_paths, "mp_before_extend");

    INFO("SUBSTAGE = mate-pair libraries ")
    exspander_stage = PathExtendStage::MPStage;
    all_libs.clear();
    max_is_right_quantile = FindOverlapLenForStage(exspander_stage);
    PathContainer mp_paths(clone_paths);

    if (is_2015_scaffolder_enabled(sc_mode)) {
        //TODO: constants
        for (auto cur_length = min_unique_length; cur_length > 500; cur_length -= 500) {
            ScaffoldingUniqueEdgeStorage current_unique_storage;
            ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp, cur_length, unique_variaton);
            unique_edge_analyzer.FillUniqueEdgeStorage(current_unique_storage);
            all_libs = MakeAllExtenders(exspander_stage, gp, clone_map, pset, current_unique_storage, clone_paths);
            shared_ptr<CompositeExtender> mp_main_pe = make_shared<CompositeExtender>(gp.g, clone_map, all_libs,
                                                                                      max_is_right_quantile,
                                                                                      main_unique_storage,
                                                                                      cfg::get().pe_params.param_set.extension_options.max_repeat_length);
            INFO("Growing paths using mate-pairs unique length " << cur_length);
            mp_paths = resolver.extendSeeds(mp_paths, *mp_main_pe);
            DebugOutputPaths(gp, output_dir, mp_paths, "mp_before_overlap_" + std::to_string(cur_length));
        }
    } else {
        all_libs = MakeAllExtenders(exspander_stage, gp, clone_map, pset, main_unique_storage, clone_paths);
        shared_ptr<CompositeExtender> mp_main_pe = make_shared<CompositeExtender>(gp.g, clone_map, all_libs,
                                                                                  max_is_right_quantile,
                                                                                  main_unique_storage,
                                                                                  cfg::get().pe_params.param_set.extension_options.max_repeat_length);
        INFO("Growing paths using mate-pairs");
        mp_paths = resolver.extendSeeds(clone_paths, *mp_main_pe);

        DebugOutputPaths(gp, output_dir, mp_paths, "mp_before_overlap");
        FinalizePaths(mp_paths, clone_map, max_is_right_quantile, max_is_right_quantile, true);
    }
    DebugOutputPaths(gp, output_dir, mp_paths, "mp_final_paths");
    DEBUG("Paths are grown with mate-pairs");

//MP end

//pe again
    INFO("SUBSTAGE = polishing paths")
    exspander_stage = PathExtendStage::FinalizingPEStage;
    all_libs.clear();
    all_libs = MakeAllExtenders(exspander_stage, gp, clone_map, pset, main_unique_storage);
    max_is_right_quantile = FindOverlapLenForStage(exspander_stage);
    shared_ptr<CompositeExtender> last_extender = make_shared<CompositeExtender>(gp.g, clone_map, all_libs,
                                                                                 max_is_right_quantile, main_unique_storage,
                                                                                 cfg::get().pe_params.param_set.extension_options.max_repeat_length);

    auto last_paths = resolver.extendSeeds(mp_paths, *last_extender);
    DebugOutputPaths(gp, output_dir, last_paths, "mp2_before_overlap");

    exspander_stage = PathExtendStage::FinalPolishing;
    all_libs = MakeAllExtenders(exspander_stage, gp, clone_map, pset, main_unique_storage);
    last_extender = make_shared<CompositeExtender>(gp.g, clone_map, all_libs,
                                            max_is_right_quantile, main_unique_storage,
                                            cfg::get().pe_params.param_set.extension_options.max_repeat_length);
    if (!is_2015_scaffolder_enabled(sc_mode)) {
        FinalizePaths(last_paths, clone_map, min_edge_len, max_is_right_quantile);
        DebugOutputPaths(gp, output_dir, last_paths, "mp2_before_traverse");
    }

    TraverseLoops(last_paths, clone_map, last_extender);
    FinalizePaths(last_paths, clone_map, min_edge_len, max_is_right_quantile);

//result
    if (broken_contigs.is_initialized()) {
        OutputBrokenScaffolds(last_paths, (int) gp.g.k(), writer, output_dir + broken_contigs.get());
    }
    debruijn_graph::GenomeConsistenceChecker genome_checker (gp, main_unique_storage, 1000, 0.2);
    DebugOutputPaths(gp, output_dir, last_paths, "mp2_final_paths");
    writer.OutputPaths(last_paths, output_dir + contigs_name);
    if (gp.genome.size() > 0)
        CountMisassembliesWithReference(genome_checker, last_paths);
    //FinalizeUniquenessPaths();

//TODO:: destructor?
    last_paths.DeleteAllPaths();
    seeds.DeleteAllPaths();
    mp_paths.DeleteAllPaths();
    clone_paths.DeleteAllPaths();

    INFO("ExSPAnder repeat resolving tool finished");
}

} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
