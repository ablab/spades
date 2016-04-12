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

struct PathExtendParamsContainer {

    PathExtendParamsContainer(const pe_config::MainPEParamsT& pe_cfg_,
                              const std::string& output_dir_,
                              const std::string& contigs_name_,
                              const std::string& scf_name_,
                              config::pipeline_type mode_,
                              bool uneven_depth_,
                              bool avoid_rc_connections_,
                              bool use_scaffolder_,
                              bool output_broken_scaffolds_ = true):
        pe_cfg(pe_cfg_),
        pset(pe_cfg_.param_set),
        output_dir(output_dir_),
        etc_dir(output_dir + pe_cfg_.etc_dir + "/"),
        contigs_name(scf_name_),
        broken_contigs(contigs_name_),
        mode(mode_),
        uneven_depth(uneven_depth_),
        avoid_rc_connections(avoid_rc_connections_),
        use_scaffolder(use_scaffolder_),
        traverse_loops(true),
        output_broken_scaffolds(output_broken_scaffolds_)
    {
        if (!(use_scaffolder && pset.scaffolder_options.enabled)) {
            contigs_name = contigs_name_;
            traverse_loops = false;
            output_broken_scaffolds = false;
        }
    }

    const pe_config::MainPEParamsT& pe_cfg;
    const pe_config::ParamSetT& pset;

    std::string output_dir;
    std::string etc_dir;

    std::string contigs_name;
    std::string broken_contigs;

    config::pipeline_type mode;
    bool uneven_depth;

    bool avoid_rc_connections;
    bool use_scaffolder;
    bool traverse_loops;
    bool output_broken_scaffolds;
};

inline void DebugOutputPaths(const conj_graph_pack& gp,
                             const PathExtendParamsContainer& params,
                             const PathContainer& paths,
                             const string& name) {
    PathInfoWriter path_writer;
    PathVisualizer visualizer;

    DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(gp.g);
    DefaultContigConstructor<ConjugateDeBruijnGraph> constructor(gp.g, corrector);
    ContigWriter writer(gp.g, constructor, gp.components, params.mode == config::pipeline_type::plasmid);

    if (!params.pe_cfg.debug_output) {
        return;
    }
    writer.OutputPaths(paths, params.etc_dir + name);
    if (params.pe_cfg.output.write_paths) {
        path_writer.WritePaths(paths, params.etc_dir + name + ".dat");
    }
    if (params.pe_cfg.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp, params.etc_dir + name + ".dot", name, paths);
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


inline void OutputBrokenScaffolds(PathContainer& paths,
                                  const PathExtendParamsContainer& params,
                                  int k,
                                  const ContigWriter& writer,
                                  const std::string& filename) {
    if (!params.pset.scaffolder_options.enabled
        || !params.use_scaffolder
        || params.pe_cfg.obs == obs_none) {
        return;
    }

    int min_gap = params.pe_cfg.obs == obs_break_all ? k / 2 : k;

    ScaffoldBreaker breaker(min_gap, paths);
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

bool HasOnlyMPLibs(const config::dataset& dataset_info) {
    for (const auto& lib : dataset_info.reads) {
        if (!((lib.type() == io::LibraryType::MatePairs || lib.type() == io::LibraryType::HQMatePairs) &&
              lib.data().mean_insert_size > 0.0)) {
            return false;
        }
    }
    return true;
}

bool UseCoverageResolverForSingleReads(const config::dataset& dataset_info,
                                       const io::LibraryType& type) {
    return HasOnlyMPLibs(dataset_info) && (type == io::LibraryType::HQMatePairs);
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

inline void FinalizePaths(const PathExtendParamsContainer& params,
                          PathContainer& paths,
                          const Graph& g,
                          GraphCoverageMap& cover_map,
                          size_t min_edge_len,
                          size_t max_path_diff,
                          bool mate_pairs = false) {
    PathExtendResolver resolver(cover_map.graph());

    if (params.pset.remove_overlaps) {
        resolver.removeOverlaps(paths, cover_map, min_edge_len, max_path_diff,
                                params.pset.cut_all_overlaps,
                                (params.mode == config::pipeline_type::moleculo));
    }
    else {
        resolver.removeEqualPaths(paths, cover_map, min_edge_len);
    }
    if (mate_pairs) {
        resolver.RemoveMatePairEnds(paths, min_edge_len);
    }
    if (params.avoid_rc_connections) {
        paths.FilterInterstandBulges();
    }
    paths.FilterEmptyPaths();
    if (!mate_pairs) {
        resolver.addUncoveredEdges(paths, cover_map);
    }
    if (params.pset.path_filtration.enabled) {
        LengthPathFilter length_filter(g, params.pset.path_filtration.min_length);
        length_filter.filter(paths);
        IsolatedPathFilter low_coverage_filter(g, params.pset.path_filtration.min_length_for_low_covered, params.pset.path_filtration.min_coverage);
        IsolatedPathFilter short_isolated_filter(g, params.pset.path_filtration.isolated_min_length);
        low_coverage_filter.filter(paths);
        short_isolated_filter.filter(paths);
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

inline bool IsForSingleReadExtender(const io::SequencingLibrary<config::DataSetData> &lib) {
    io::LibraryType lt = lib.type();
    return (lib.data().single_reads_mapped ||
            lt == io::LibraryType::PacBioReads ||
            lt == io::LibraryType::SangerReads ||
            lt == io::LibraryType::NanoporeReads ||
            lib.is_contig_lib());
}

inline bool IsForPEExtender(const io::SequencingLibrary<config::DataSetData> &lib) {
    return (lib.type() == io::LibraryType::PairedEnd &&
            lib.data().mean_insert_size > 0.0);
}

inline bool IsForShortLoopExtender(const io::SequencingLibrary<config::DataSetData> &lib) {
    return (lib.type() == io::LibraryType::PairedEnd &&
            lib.data().mean_insert_size > 0.0);
}

inline bool IsForScaffoldingExtender(const io::SequencingLibrary<config::DataSetData> &lib) {
    return (lib.type() == io::LibraryType::PairedEnd &&
            lib.data().mean_insert_size > 0.0);
}

inline bool IsForMPExtender(const io::SequencingLibrary<config::DataSetData> &lib) {
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
inline shared_ptr<PairedInfoLibrary> MakeNewLib(const config::dataset::Library& lib,
                                                const conj_graph_pack::graph_t& g,
                                                const Index& paired_index) {
    size_t read_length = lib.data().read_length;
    size_t is = (size_t) lib.data().mean_insert_size;
    int is_min = (int) lib.data().insert_size_left_quantile;
    int is_max = (int) lib.data().insert_size_right_quantile;
    int var = (int) lib.data().insert_size_deviation;
    bool is_mp = lib.type() == io::LibraryType::MatePairs ||  lib.type() == io::LibraryType::HQMatePairs ;
    return make_shared< PairedInfoLibraryWithIndex<decltype(paired_index)> >(g.k(), g, read_length,
                                                                                    is, is_min > 0.0 ? size_t(is_min) : 0, is_max > 0.0 ? size_t(is_max) : 0,
                                                                                    size_t(var),
                                                                                    paired_index, is_mp,
                                                                                    lib.data().insert_size_distribution);
}

pe_config::LongReads GetLongReadsConfig(const PathExtendParamsContainer& params,
                                        const io::LibraryType& type) {
    if (io::SequencingLibraryBase::is_long_read_lib(type)) {
        return params.pe_cfg.long_reads.pacbio_reads;
    } else if (type == io::LibraryType::PathExtendContigs){
        return params.pe_cfg.long_reads.meta_contigs;
    } else if (io::SequencingLibraryBase::is_contig_lib(type)) {
        return params.pe_cfg.long_reads.contigs;
    }
    return params.pe_cfg.long_reads.single_reads;
}


inline shared_ptr<ExtensionChooser> MakeLongReadsExtensionChooser(const config::dataset::Library& lib,
                                                                  size_t lib_index,
                                                                  const PathExtendParamsContainer& params,
                                                                  const conj_graph_pack& gp) {
    PathContainer paths;
    AddPathsToContainer(gp, gp.single_long_reads[lib_index].GetAllPaths(), 1, paths);

    auto long_reads_config = GetLongReadsConfig(params, lib.type());
    return make_shared<LongReadsExtensionChooser>(gp.g, paths, long_reads_config.filtering,
                                                  long_reads_config.weight_priority,
                                                  long_reads_config.unique_edge_priority,
                                                  long_reads_config.min_significant_overlap,
                                                  params.pset.extension_options.max_repeat_length,
                                                  params.uneven_depth);
}


inline shared_ptr<SimpleExtender> MakeLongReadsExtender(const config::dataset& dataset_info,
                                                        size_t lib_index,
                                                        const PathExtendParamsContainer& params,
                                                        const conj_graph_pack& gp,
                                                        const GraphCoverageMap& cov_map) {
    const auto& lib = dataset_info.reads[lib_index];
    size_t resolvable_repeat_length_bound = 10000ul;
    if (!lib.is_contig_lib()) {
        resolvable_repeat_length_bound = std::max(resolvable_repeat_length_bound, lib.data().read_length);
    }
    INFO("resolvable_repeat_length_bound set to " << resolvable_repeat_length_bound);


    auto long_read_ec = MakeLongReadsExtensionChooser(lib, lib_index, params, gp);
    return make_shared<SimpleExtender>(gp, cov_map,
                                       long_read_ec,
                                       resolvable_repeat_length_bound,
                                       params.pset.loop_removal.max_loops,
                                       true, /* investigate short loops */
                                       UseCoverageResolverForSingleReads(dataset_info, lib.type()));
}

inline shared_ptr<SimpleExtender> MakeLongEdgePEExtender(const config::dataset& dataset_info,
                                                         size_t lib_index,
                                                         const PathExtendParamsContainer& params,
                                                         const conj_graph_pack& gp,
                                                         const GraphCoverageMap& cov_map,
                                                         bool investigate_loops) {

    const auto& lib = dataset_info.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(lib, gp.g, gp.clustered_indices[lib_index]);
    SetSingleThresholdForLib(paired_lib, params.pset, lib.data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<WeightCounter> wc =
        make_shared<PathCoverWeightCounter>(gp.g, paired_lib, params.pset.normalize_weight);
    shared_ptr<ExtensionChooser> extension =
        make_shared<LongEdgeExtensionChooser>(gp.g, wc,
                                              GetWeightThreshold(paired_lib, params.pset),
                                              GetPriorityCoeff(paired_lib, params.pset));

    return make_shared<SimpleExtender>(gp, cov_map,
                                       extension,
                                       paired_lib->GetISMax(),
                                       params.pset.loop_removal.max_loops,
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/);
}

inline shared_ptr<SimpleExtensionChooser> MakeMetaExtensionChooser(shared_ptr<PairedInfoLibrary> lib,
                                                                   const PathExtendParamsContainer& params,
                                                                   const conj_graph_pack& gp,
                                                                   size_t read_length) {
    VERIFY(params.mode == config::pipeline_type::meta);
    VERIFY(!lib->IsMp());
    shared_ptr<WeightCounter> wc = make_shared<MetagenomicWeightCounter>(gp.g,
                                                                         lib,
                                                                         read_length, //read_length
                                                                         0.3, //normalized_threshold
                                                                         3, //raw_threshold
                                                                         0 /*estimation_edge_length*/ );
    return make_shared<SimpleExtensionChooser>(gp.g, wc,
                                               params.pset.extension_options.weight_threshold,
                                               params.pset.extension_options.priority_coeff);
}

inline shared_ptr<SimpleExtender> MakeMetaExtender(const config::dataset& dataset_info,
                                                   size_t lib_index,
                                                   const PathExtendParamsContainer& params,
                                                   const conj_graph_pack& gp,
                                                   const GraphCoverageMap& cov_map,
                                                   bool investigate_loops) {

    const auto& lib = dataset_info.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(lib, gp.g, gp.clustered_indices[lib_index]);

    return make_shared<SimpleExtender>(gp, cov_map,
                                       MakeMetaExtensionChooser(paired_lib, params, gp, dataset_info.RL()),
                                       paired_lib->GetISMax(),
                                       params.pset.loop_removal.max_loops,
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/);
}

inline shared_ptr<SimpleExtender> MakePEExtender(const config::dataset& dataset_info,
                                                 size_t lib_index,
                                                 const PathExtendParamsContainer& params,
                                                 const conj_graph_pack& gp,
                                                 const GraphCoverageMap& cov_map,
                                                 bool investigate_loops) {

    const auto& lib = dataset_info.reads[lib_index];
    shared_ptr<PairedInfoLibrary>  paired_lib = MakeNewLib(lib, gp.g, gp.clustered_indices[lib_index]);
    SetSingleThresholdForLib(paired_lib, params.pset, lib.data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<WeightCounter> wc = make_shared<PathCoverWeightCounter>(gp.g, paired_lib, params.pset.normalize_weight);
    auto extension = make_shared<SimpleExtensionChooser>(gp.g, wc,
                                                         GetWeightThreshold(paired_lib, params.pset),
                                                         GetPriorityCoeff(paired_lib, params.pset));

    return make_shared<SimpleExtender>(gp, cov_map,
                                       extension,
                                       paired_lib->GetISMax(),
                                       params.pset.loop_removal.max_loops,
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/);
}


inline shared_ptr<PathExtender> MakeScaffoldingExtender(const config::dataset& dataset_info,
                                                        size_t lib_index,
                                                        const PathExtendParamsContainer& params,
                                                        const conj_graph_pack& gp,
                                                        const GraphCoverageMap& cov_map) {
    const auto& lib = dataset_info.reads[lib_index];
    const auto& pset = params.pset;
    shared_ptr<PairedInfoLibrary>  paired_lib = MakeNewLib(lib, gp.g, gp.scaffolding_indices[lib_index]);

    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(gp.g, paired_lib);

    auto scaff_chooser = std::make_shared<ScaffoldingExtensionChooser>(gp.g, counter,
                                                                       pset.scaffolder_options.cl_threshold,
                                                                       pset.scaffolder_options.var_coeff);

    vector<shared_ptr<GapJoiner>> joiners;
    if (params.pset.scaffolder_options.use_la_gap_joiner)
        joiners.push_back(std::make_shared<LAGapJoiner>(gp.g, pset.scaffolder_options.min_overlap_length,
                                                        pset.scaffolder_options.flank_multiplication_coefficient,
                                                        pset.scaffolder_options.flank_addition_coefficient));


    joiners.push_back(std::make_shared<HammingGapJoiner>(gp.g,
                                                         pset.scaffolder_options.min_gap_score,
                                                         pset.scaffolder_options.short_overlap,
                                                         (int) pset.scaffolder_options.basic_overlap_coeff * dataset_info.RL()));

    auto composite_gap_joiner = std::make_shared<CompositeGapJoiner>(gp.g, 
                                                joiners, 
                                                size_t(pset.scaffolder_options.max_can_overlap * (double) gp.g.k()), /* may overlap threshold */
                                                int(math::round((double) gp.g.k() - pset.scaffolder_options.var_coeff * (double) paired_lib->GetIsVar())),  /* must overlap threshold */
                                                pset.scaffolder_options.artificial_gap);

    return make_shared<ScaffoldingPathExtender>(gp, cov_map, scaff_chooser,
                                                composite_gap_joiner,
                                                paired_lib->GetISMax(),
                                                pset.loop_removal.max_loops,
                                                false, /* investigate short loops */
                                                params.avoid_rc_connections);
}


inline shared_ptr<PathExtender> MakeRNAScaffoldingExtender(const config::dataset& dataset_info,
                                                            size_t lib_index,
                                                            const PathExtendParamsContainer& params,
                                                            const conj_graph_pack& gp,
                                                            const GraphCoverageMap& cov_map) {

    const auto& lib = dataset_info.reads[lib_index];
    const auto& pset = params.pset;
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(lib, gp.g, gp.paired_indices[lib_index]);

    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(gp.g, paired_lib);

    auto scaff_chooser = std::make_shared<ScaffoldingExtensionChooser>(gp.g, counter, pset.scaffolder_options.cutoff, pset.scaffolder_options.var_coeff);
    auto scaff_chooser2 = std::make_shared<ScaffoldingExtensionChooser>(gp.g, counter, pset.scaffolder_options.hard_cutoff, pset.scaffolder_options.var_coeff);

    vector<shared_ptr<GapJoiner>> joiners;
    if (params.pset.scaffolder_options.use_la_gap_joiner)
        joiners.push_back(std::make_shared<LAGapJoiner>(gp.g, pset.scaffolder_options.min_overlap_length,
                                                        pset.scaffolder_options.flank_multiplication_coefficient,
                                                        pset.scaffolder_options.flank_addition_coefficient));


    joiners.push_back(std::make_shared<HammingGapJoiner>(gp.g,
                                                         pset.scaffolder_options.min_gap_score,
                                                         pset.scaffolder_options.short_overlap,
                                                         (int) pset.scaffolder_options.basic_overlap_coeff * dataset_info.RL()));

    auto composite_gap_joiner = std::make_shared<CompositeGapJoiner>(gp.g,
                                                                     joiners,
                                                                     size_t(pset.scaffolder_options.max_can_overlap * (double) gp.g.k()), /* may overlap threshold */
                                                                     int(math::round((double) gp.g.k() - pset.scaffolder_options.var_coeff * (double) paired_lib->GetIsVar())),  /* must overlap threshold */
                                                                     pset.scaffolder_options.artificial_gap);

    VERIFY(pset.scaffolder_options.min_overlap_for_rna_scaffolding.is_initialized());
    return make_shared<RNAScaffoldingPathExtender>(gp, cov_map,
                                                   scaff_chooser,
                                                   scaff_chooser2,
                                                   composite_gap_joiner,
                                                   paired_lib->GetISMax(),
                                                   pset.loop_removal.max_loops,
                                                   false  /* investigate short loops */,
                                                   *pset.scaffolder_options.min_overlap_for_rna_scaffolding);
}


inline shared_ptr<PathExtender> MakeScaffolding2015Extender(const config::dataset& dataset_info,
                                                            size_t lib_index,
                                                            const PathExtendParamsContainer& params,
                                                            const conj_graph_pack& gp,
                                                            const GraphCoverageMap& cov_map,
                                                            const ScaffoldingUniqueEdgeStorage& storage) {

    const auto& lib = dataset_info.reads[lib_index];
    const auto& pset = params.pset;
    shared_ptr<PairedInfoLibrary> paired_lib;
    INFO("Creating Scaffolding 2015 extender for lib #" << lib_index);

    //TODO:: temporary solution
    if (gp.paired_indices[lib_index].size() > gp.clustered_indices[lib_index].size()) {
        INFO("Paired unclustered indices not empty, using them");
        paired_lib = MakeNewLib(lib, gp.g, gp.paired_indices[lib_index]);
    } else if (gp.clustered_indices[lib_index].size() != 0 ) {
        INFO("clustered indices not empty, using them");
        paired_lib = MakeNewLib(lib, gp.g, gp.clustered_indices[lib_index]);
    } else {
        ERROR("All paired indices are empty!");
    }

    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(gp.g, paired_lib);
//TODO::was copypasted from MakeScaffoldingExtender, refactor 2015 extension chhoser
    DEBUG("creating extchooser");

    auto scaff_chooser = std::make_shared<ExtensionChooser2015>(gp.g,
                                                                counter,
                                                                lib_index,
                                                                storage,
                                                                pset.scaffolder_options.cl_threshold,
                                                                pset.scaffolder_options.var_coeff,
                                                                pset.scaffolding2015.relative_weight_cutoff);

    auto gap_joiner = std::make_shared<HammingGapJoiner>(gp.g, pset.scaffolder_options.min_gap_score,
                                                         pset.scaffolder_options.short_overlap,
                                                         (int) pset.scaffolder_options.basic_overlap_coeff * dataset_info.RL());

    return make_shared<ScaffoldingPathExtender>(gp, cov_map,
                                                scaff_chooser,
                                                gap_joiner,
                                                paired_lib->GetISMax(),
                                                pset.loop_removal.max_loops,
                                                false, /* investigate short loops */
                                                params.avoid_rc_connections,
                                                false /* jump only from tips */);
}


inline shared_ptr<SimpleExtender> MakeMPExtender(const config::dataset& dataset_info,
                                                 size_t lib_index,
                                                 const PathExtendParamsContainer& params,
                                                 const conj_graph_pack& gp,
                                                 const GraphCoverageMap& cov_map,
                                                 const PathContainer& paths) {

    const auto& lib = dataset_info.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(lib, gp.g, gp.paired_indices[lib_index]);

    SetSingleThresholdForLib(paired_lib, params.pset, lib.data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    size_t max_number_of_paths_to_search = GetNumberMPPaths(gp.g);
    DEBUG("max number of mp paths " << max_number_of_paths_to_search);

    shared_ptr<MatePairExtensionChooser> chooser =
        make_shared<MatePairExtensionChooser>(gp.g,
                                              paired_lib,
                                              paths,
                                              max_number_of_paths_to_search,
                                              params.uneven_depth);

    return make_shared<SimpleExtender>(gp, cov_map,
                                       chooser,
                                       paired_lib->GetISMax(),
                                       params.pset.loop_removal.mp_max_loops,
                                       true, /* investigate short loops */
                                       false /*use short loop coverage resolver*/);
}


inline shared_ptr<SimpleExtender> MakeCoordCoverageExtender(const config::dataset& dataset_info,
                                                            size_t lib_index,
                                                            const PathExtendParamsContainer& params,
                                                            const conj_graph_pack& gp,
                                                            const GraphCoverageMap& cov_map) {

    const auto& lib = dataset_info.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(lib, gp.g, gp.clustered_indices[lib_index]);

    CoverageAwareIdealInfoProvider provider(gp.g, paired_lib, -1ul, 0);
    auto coord_chooser = make_shared<CoordinatedCoverageExtensionChooser>(gp.g, provider,
                                                                          params.pset.coordinated_coverage.max_edge_length_in_repeat,
                                                                          params.pset.coordinated_coverage.delta,
                                                                          params.pset.coordinated_coverage.min_path_len);
    auto chooser = make_shared<JointExtensionChooser>(gp.g, MakeMetaExtensionChooser(paired_lib, params, gp, dataset_info.RL()), coord_chooser);

    return make_shared<SimpleExtender>(gp, cov_map, chooser,
                                       -1ul /* insert size */,
                                       params.pset.loop_removal.mp_max_loops,
                                       true, /* investigate short loops */
                                       false /*use short loop coverage resolver*/);
}


inline shared_ptr<SimpleExtender> MakeRNAExtender(const config::dataset& dataset_info,
                                                  size_t lib_index,
                                                  const PathExtendParamsContainer& params,
                                                  const conj_graph_pack& gp,
                                                  const GraphCoverageMap& cov_map,
                                                  bool investigate_loops) {

    const auto& lib = dataset_info.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(lib, gp.g, gp.clustered_indices[lib_index]);
    SetSingleThresholdForLib(paired_lib, params.pset, lib.data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<WeightCounter> wc = make_shared<PathCoverWeightCounter>(gp.g, paired_lib, params.pset.normalize_weight);
    shared_ptr<RNAExtensionChooser> extension =
        make_shared<RNAExtensionChooser>(gp.g, wc,
                                         GetWeightThreshold(paired_lib, params.pset),
                                         GetPriorityCoeff(paired_lib, params.pset));

    return make_shared<MultiExtender>(gp, cov_map,
                                      extension,
                                      paired_lib->GetISMax(),
                                      params.pset.loop_removal.max_loops,
                                      investigate_loops,
                                      false /*use short loop coverage resolver*/);
}


inline shared_ptr<SimpleExtender> MakeRNALongReadsExtender(const config::dataset& dataset_info,
                                                           size_t lib_index,
                                                           const PathExtendParamsContainer& params,
                                                           const conj_graph_pack& gp,
                                                           const GraphCoverageMap& cov_map) {

    VERIFY_MSG(false, "Long reads rna extender is not implemented yet")

    const auto& lib = dataset_info.reads[lib_index];
    size_t resolvable_repeat_length_bound = 10000ul;
    if (!lib.is_contig_lib()) {
        resolvable_repeat_length_bound = std::max(resolvable_repeat_length_bound, lib.data().read_length);
    }
    INFO("resolvable_repeat_length_bound set to " << resolvable_repeat_length_bound);

    auto long_reads_ec = MakeLongReadsExtensionChooser(lib, lib_index, params, gp);

    return make_shared<SimpleExtender>(gp, cov_map,
                                       long_reads_ec,
                                       resolvable_repeat_length_bound,
                                       params.pset.loop_removal.max_loops,
                                       true, /* investigate short loops */
                                       UseCoverageResolverForSingleReads(dataset_info, lib.type()));
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

inline vector<shared_ptr<PathExtender> > MakeAllExtenders(PathExtendStage stage,
                                                          const config::dataset& dataset_info,
                                                          const PathExtendParamsContainer& params,
                                                          const conj_graph_pack& gp,
                                                          const GraphCoverageMap& cov_map,
                                                          const ScaffoldingUniqueEdgeStorage& storage,
                                                          const PathContainer& paths_for_mp = PathContainer()) {

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

    const auto& pset = params.pset;

    for (io::LibraryType lt : io::LibraryPriotity) {
        for (size_t lib_index = 0; lib_index < dataset_info.reads.lib_count(); ++lib_index) {
            const auto& lib = dataset_info.reads[lib_index];
            if (lib.type() != lt)
                continue;

            //TODO: scaff2015 does not need any single read libs?
            if (IsForSingleReadExtender(lib) && pset.sm != sm_2015) {
                result.push_back(MakeLongReadsExtender(dataset_info, lib_index, params, gp, cov_map));
                ++single_read_libs;
            }
            if (IsForPEExtender(lib)) {
                ++pe_libs;
                if (IsPEStage(stage) && IsOldPEEnabled(pset.sm)) {
                    if (params.mode == config::pipeline_type::meta)
                        //TODO proper configuration via config
                        pes.push_back(MakeMetaExtender(dataset_info, lib_index, params, gp, cov_map, false));
                    else if (params.mode == config::pipeline_type::moleculo)
                        pes.push_back(MakeLongEdgePEExtender(dataset_info, lib_index, params, gp, cov_map, false));
                    else if (pset.multi_path_extend  && !IsPolishingStage(stage))
                        pes.push_back(MakeRNAExtender(dataset_info, lib_index, params, gp, cov_map, false));
                    else
                        pes.push_back(MakePEExtender(dataset_info, lib_index, params, gp, cov_map, false));
                }
                else if (pset.sm == sm_2015) {
                    pes2015.push_back(MakeScaffolding2015Extender(dataset_info, lib_index, params, gp, cov_map, storage));
                }
            }
            //FIXME logic is very cryptic!
            if (IsForShortLoopExtender(lib) && IsOldPEEnabled(pset.sm)) {
                if (params.mode == config::pipeline_type::meta)
                    pes.push_back(MakeMetaExtender(dataset_info, lib_index, params, gp, cov_map, true));
                else if (pset.multi_path_extend && !IsPolishingStage(stage))
                    pes.push_back(MakeRNAExtender(dataset_info, lib_index, params, gp, cov_map, true));
                else
                    pe_loops.push_back(MakePEExtender(dataset_info, lib_index, params, gp, cov_map, true));
            }
            if (IsForScaffoldingExtender(lib) && params.use_scaffolder && pset.scaffolder_options.enabled) {
                ++scf_pe_libs;
                if (params.mode == config::pipeline_type::rna) {
                    pe_scafs.push_back(MakeRNAScaffoldingExtender(dataset_info, lib_index, params, gp, cov_map));
                }
                else {
                    switch (pset.sm) {
                        case sm_old: {
                            pe_scafs.push_back(MakeScaffoldingExtender(dataset_info, lib_index, params, gp, cov_map));
                            break;
                        }
                        case sm_old_pe_2015: {
                            pe_scafs.push_back(MakeScaffolding2015Extender(dataset_info, lib_index, params, gp, cov_map, storage));
                            break;
                        }
                        case sm_combined: {
                            pe_scafs.push_back(MakeScaffoldingExtender(dataset_info, lib_index, params, gp, cov_map));
                            pe_scafs.push_back(MakeScaffolding2015Extender(dataset_info, lib_index, params, gp, cov_map, storage));
                            break;
                        }
                        default:
                            break;
                    }
                }
            }
            if (IsForMPExtender(lib) && IsMPStage(stage)) {
                ++mp_libs;
                switch (pset.sm) {
                    case sm_old: {
                        mps.push_back(MakeMPExtender(dataset_info, lib_index, params, gp, cov_map, paths_for_mp));
                        break;
                    }
                    case sm_old_pe_2015: {
                        mps.push_back(MakeScaffolding2015Extender(dataset_info, lib_index, params, gp, cov_map, storage));
                        break;
                    }
                    case sm_2015: {
                        mps.push_back(MakeScaffolding2015Extender(dataset_info, lib_index, params, gp, cov_map, storage));
                        break;
                    }
                    case sm_combined: {
                        mps.push_back(MakeMPExtender(dataset_info, lib_index, params, gp, cov_map, paths_for_mp));
                        mps.push_back(MakeScaffolding2015Extender(dataset_info, lib_index, params, gp, cov_map, storage));
                        break;
                    }
                    default:
                        break;
                }
            }
        }

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
    INFO("Scaffolder is " << (pset.scaffolder_options.enabled ? "on" : "off"));

    if (pset.use_coordinated_coverage) {
        INFO("Using additional coordinated coverage extender");
        result.push_back(MakeCoordCoverageExtender(dataset_info, 0 /* lib index */, params, gp, cov_map));
    }

    PrintExtenders(result);
    return result;
}

inline shared_ptr<scaffold_graph::ScaffoldGraph> ConstructScaffoldGraph(const config::dataset& dataset_info,
                                                                        const pe_config::ParamSetT::ScaffoldGraphParamsT& params,
                                                                        const conj_graph_pack& gp,
                                                                        const ScaffoldingUniqueEdgeStorage& edge_storage) {
    using namespace scaffold_graph;
    vector<shared_ptr<ConnectionCondition>> conditions;

    INFO("Constructing connections");
    LengthEdgeCondition edge_condition(gp.g, edge_storage.GetMinLength());

    for (size_t lib_index = 0; lib_index < dataset_info.reads.lib_count(); ++lib_index) {
        const auto& lib = dataset_info.reads[lib_index];
        if (lib.is_paired()) {
            shared_ptr<PairedInfoLibrary> paired_lib;
            if (IsForMPExtender(lib))
                paired_lib = MakeNewLib(lib, gp.g, gp.paired_indices[lib_index]);
            else if (IsForPEExtender(lib))
                paired_lib = MakeNewLib(lib, gp.g, gp.clustered_indices[lib_index]);
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
                               const set<EdgeId>& main_edge_set,
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


inline size_t FindOverlapLenForStage(PathExtendStage stage, const config::dataset& dataset_info) {
    size_t res = 0;
    for (const auto& lib : dataset_info.reads) {
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

inline bool MPLibsExist(const config::dataset& dataset_info) {
    for (const auto& lib : dataset_info.reads)
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
                                                          const config::dataset& dataset_info,
                                                          size_t& min_unique_length,
                                                          double& unique_variation,
                                                          bool autodetect) {

    ScaffoldingUniqueEdgeStorage main_unique_storage;
    //Setting scaffolding2015 parameters
    if (autodetect) {
        INFO("Autodetecting unique edge set parameters...");
        bool pe_found = false;
        //TODO constants
        size_t min_MP_IS = 10000;
        for (size_t i = 0; i < dataset_info.reads.lib_count(); ++i) {

            if (IsForPEExtender(dataset_info.reads[i])) {
                pe_found = true;
            }
            if (IsForMPExtender(dataset_info.reads[i])) {
                min_MP_IS = min(min_MP_IS, (size_t) dataset_info.reads[i].data().mean_insert_size);
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


inline void ResolveRepeatsPe(const config::dataset& dataset_info,
                             const PathExtendParamsContainer& params,
                             conj_graph_pack& gp) {

    INFO("ExSPAnder repeat resolving tool started");
    const pe_config::ParamSetT &pset = params.pset;

    ScaffoldingUniqueEdgeStorage main_unique_storage;
    auto sc_mode = pset.sm;
    auto min_unique_length = pset.scaffolding2015.min_unique_length;
    auto unique_variaton = pset.scaffolding2015.unique_coverage_variation;
    bool detect_repeats_online = !(IsScaffolder2015Enabled(sc_mode) || params.mode == config::pipeline_type::meta);

    //Fill the storage to enable unique edge check
    if (IsScaffolder2015Enabled(sc_mode)) {
        main_unique_storage = FillUniqueEdgeStorage(gp, dataset_info,
                                                    min_unique_length,
                                                    unique_variaton,
                                                    pset.scaffolding2015.autodetect);
    }

    make_dir(params.output_dir);
    make_dir(params.etc_dir);


    //Scaffold graph
    shared_ptr<scaffold_graph::ScaffoldGraph> scaffoldGraph;
    if (pset.scaffold_graph_params.construct) {
        scaffoldGraph = ConstructScaffoldGraph(dataset_info, params.pset.scaffold_graph_params, gp, main_unique_storage);
        if (pset.scaffold_graph_params.output) {
            PrintScaffoldGraph(scaffoldGraph, main_unique_storage.GetSet(), params.etc_dir + "scaffold_graph");
        }
    }


    DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(gp.g);
    DefaultContigConstructor<ConjugateDeBruijnGraph> constructor(gp.g, corrector);
    ContigWriter writer(gp.g, constructor, gp.components, params.mode == config::pipeline_type::plasmid);


//make pe + long reads extenders
    GraphCoverageMap cover_map(gp.g);
    INFO("SUBSTAGE = paired-end libraries")
    PathExtendStage exspander_stage = PathExtendStage::PEStage;
    vector<shared_ptr<PathExtender> > all_libs =
        MakeAllExtenders(exspander_stage, dataset_info, params, gp, cover_map, main_unique_storage);

    //Parameters are subject to change
    size_t max_is_right_quantile = max(FindOverlapLenForStage(exspander_stage, dataset_info), gp.g.k() + 100);
    size_t min_edge_len = 100;
    size_t max_edge_diff_pe = /*cfg::get().mode == config::pipeline_type::rna ? 0 :*/ max_is_right_quantile;

    shared_ptr<CompositeExtender> mainPE = make_shared<CompositeExtender>(gp.g, cover_map, all_libs,
                                                                          main_unique_storage,
                                                                          max_is_right_quantile,
                                                                          pset.extension_options.max_repeat_length,
                                                                          detect_repeats_online);

//extend pe + long reads
    PathExtendResolver resolver(gp.g);
    auto seeds = resolver.makeSimpleSeeds();
    DebugOutputPaths(gp, params, seeds, "init_paths");
    seeds.SortByLength();
    INFO("Growing paths using paired-end and long single reads");
    INFO("Multi path extend is " << (cfg::get().pe_params.param_set.multi_path_extend ? "on" : "off"))
    INFO("Overlap removal is " << (cfg::get().pe_params.param_set.remove_overlaps ? "on" : "off"))
    auto paths = resolver.extendSeeds(seeds, *mainPE);
    paths.SortByLength();
    DebugOutputPaths(gp, params, paths, "pe_before_overlap");

    PathContainer clone_paths;
    GraphCoverageMap clone_map(gp.g);
    bool mp_exist = MPLibsExist(dataset_info);

    if (mp_exist) {
        ClonePathContainer(paths, clone_paths, clone_map);
    }

    exspander_stage = PathExtendStage::PEPolishing;
    all_libs = MakeAllExtenders(exspander_stage, dataset_info, params, gp, cover_map, main_unique_storage);
    mainPE = make_shared<CompositeExtender>(gp.g, cover_map, all_libs,
                                            main_unique_storage,
                                            max_is_right_quantile,
                                            pset.extension_options.max_repeat_length,
                                            detect_repeats_online);

    //We do not run overlap removal in 2015 mode
    if (!IsScaffolder2015Enabled(sc_mode))
        FinalizePaths(params, paths, gp.g, cover_map, min_edge_len, max_edge_diff_pe);
    if (params.output_broken_scaffolds) {
        OutputBrokenScaffolds(paths, params, (int) gp.g.k(), writer,
                              params.output_dir + (mp_exist ? "pe_contigs" : params.broken_contigs));
    }
    DebugOutputPaths(gp, params, paths, "pe_before_traverse");
    if (params.traverse_loops) {
        TraverseLoops(paths, cover_map, mainPE);
        FinalizePaths(params, paths, gp.g, cover_map, min_edge_len, max_edge_diff_pe);
    }
    DebugOutputPaths(gp, params, paths, (mp_exist ? "pe_final_paths" : "final_paths"));
    writer.OutputPaths(paths, params.output_dir + (mp_exist ? "pe_scaffolds" : params.contigs_name));

    cover_map.Clear();
    seeds.DeleteAllPaths();
    paths.DeleteAllPaths();
    if (!mp_exist) {
        return;
    }

//MP
    DebugOutputPaths(gp, params, clone_paths, "mp_before_extend");

    INFO("SUBSTAGE = mate-pair libraries ")
    exspander_stage = PathExtendStage::MPStage;
    all_libs.clear();
    max_is_right_quantile = FindOverlapLenForStage(exspander_stage, dataset_info);
    PathContainer mp_paths(clone_paths);

    if (IsScaffolder2015Enabled(sc_mode)) {
        //TODO: constants
        for (auto cur_length = min_unique_length; cur_length > 500; cur_length -= 500) {
            ScaffoldingUniqueEdgeStorage current_unique_storage;
            ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp, cur_length, unique_variaton);
            unique_edge_analyzer.FillUniqueEdgeStorage(current_unique_storage);
            all_libs = MakeAllExtenders(exspander_stage, dataset_info, params, gp, clone_map, current_unique_storage, clone_paths);
            shared_ptr<CompositeExtender> mp_main_pe = make_shared<CompositeExtender>(gp.g, clone_map, all_libs,
                                                                                      main_unique_storage,
                                                                                      max_is_right_quantile,
                                                                                      pset.extension_options.max_repeat_length,
                                                                                      detect_repeats_online);
            INFO("Growing paths using mate-pairs unique length " << cur_length);
            mp_paths = resolver.extendSeeds(mp_paths, *mp_main_pe);
            DebugOutputPaths(gp, params, mp_paths, "mp_before_overlap_" + std::to_string(cur_length));
        }
    } else {
        all_libs = MakeAllExtenders(exspander_stage, dataset_info, params, gp, clone_map, main_unique_storage, clone_paths);
        shared_ptr<CompositeExtender> mp_main_pe = make_shared<CompositeExtender>(gp.g, clone_map, all_libs,
                                                                                  main_unique_storage,
                                                                                  max_is_right_quantile,
                                                                                  pset.extension_options.max_repeat_length,
                                                                                  detect_repeats_online);
        INFO("Growing paths using mate-pairs");
        mp_paths = resolver.extendSeeds(clone_paths, *mp_main_pe);

        DebugOutputPaths(gp, params, mp_paths, "mp_before_overlap");
        FinalizePaths(params, mp_paths, gp.g, clone_map, max_is_right_quantile, max_is_right_quantile, true);
    }
    DebugOutputPaths(gp, params, mp_paths, "mp_final_paths");
    DEBUG("Paths are grown with mate-pairs");

//MP end

//pe again
    INFO("SUBSTAGE = polishing paths")
    exspander_stage = PathExtendStage::FinalizingPEStage;
    all_libs.clear();
    all_libs = MakeAllExtenders(exspander_stage, dataset_info, params, gp, clone_map, main_unique_storage);
    max_is_right_quantile = FindOverlapLenForStage(exspander_stage, dataset_info);
    shared_ptr<CompositeExtender> last_extender = make_shared<CompositeExtender>(gp.g, clone_map, all_libs,
                                                                                 main_unique_storage,
                                                                                 max_is_right_quantile,
                                                                                 pset.extension_options.max_repeat_length,
                                                                                 detect_repeats_online);

    auto last_paths = resolver.extendSeeds(mp_paths, *last_extender);
    DebugOutputPaths(gp, params, last_paths, "mp2_before_overlap");

    exspander_stage = PathExtendStage::FinalPolishing;
    all_libs = MakeAllExtenders(exspander_stage, dataset_info, params, gp, clone_map, main_unique_storage);
    last_extender = make_shared<CompositeExtender>(gp.g, clone_map, all_libs,
                                                   main_unique_storage,
                                                   max_is_right_quantile,
                                                   pset.extension_options.max_repeat_length,
                                                   detect_repeats_online);
    if (!IsScaffolder2015Enabled(sc_mode)) {
        FinalizePaths(params, last_paths, gp.g, clone_map, min_edge_len, max_is_right_quantile);
        DebugOutputPaths(gp, params, last_paths, "mp2_before_traverse");
    }

    TraverseLoops(last_paths, clone_map, last_extender);
    FinalizePaths(params, last_paths, gp.g, clone_map, min_edge_len, max_is_right_quantile);

//result
    if (params.output_broken_scaffolds) {
        OutputBrokenScaffolds(last_paths, params, (int) gp.g.k(), writer, params.output_dir + params.broken_contigs);
    }
    debruijn_graph::GenomeConsistenceChecker genome_checker (gp, main_unique_storage, 1000, 0.2);
    DebugOutputPaths(gp, params, last_paths, "mp2_final_paths");
    writer.OutputPaths(last_paths, params.output_dir + params.contigs_name);
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
