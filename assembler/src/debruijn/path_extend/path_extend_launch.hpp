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

#include "pe_config_struct.hpp"
#include "pe_resolver.hpp"
#include "path_extender.hpp"
#include "pe_io.hpp"
#include "path_visualizer.hpp"
#include "loop_traverser.hpp"
#include "long_read_storage.hpp"
#include "next_path_searcher.hpp"
#include "extension_chooser2015.hpp"
#include "../genome_consistance_checker.hpp"

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

inline void DebugOutputPaths(const ContigWriter& writer, const conj_graph_pack& gp,
                      const std::string& output_dir, const PathContainer& paths,
                      const string& name) {
    PathInfoWriter path_writer;
    PathVisualizer visualizer;
    string etcDir = GetEtcDir(output_dir);
    if (!cfg::get().pe_params.debug_output) {
        return;
    }
    if (cfg::get().pe_params.output.write_paths) {
        writer.writePathEdges(paths, etcDir + name + ".dat");
    }
    if (cfg::get().pe_params.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + name + ".dot", name,
                                             paths);
        path_writer.writePaths(paths, etcDir + name + ".data");
    }
}

inline double GetWeightThreshold(shared_ptr<PairedInfoLibrary> lib, const pe_config::ParamSetT& pset) {
    return lib->IsMp() ? pset.mate_pair_options.weight_threshold : pset.extension_options.weight_threshold;
}

inline double GetPriorityCoeff(shared_ptr<PairedInfoLibrary> lib, const pe_config::ParamSetT& pset) {
    return lib->IsMp() ? pset.mate_pair_options.priority_coeff : pset.extension_options.priority_coeff;
}

inline void SetSingleThresholdForLib(shared_ptr<PairedInfoLibrary> lib, const pe_config::ParamSetT &pset, double threshold) {
    if  (lib->IsMp())
        lib->SetSingleThreshold(pset.mate_pair_options.use_default_single_threshold || math::le(threshold, 0.0) ?
                                pset.mate_pair_options.single_threshold : threshold);
    else
        lib->SetSingleThreshold(pset.extension_options.use_default_single_threshold || math::le(threshold, 0.0) ?
                                pset.extension_options.single_threshold : threshold);
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
    //writer.writePaths(breaker.container(), filename + ".fasta");
    writer.WritePathsToFASTG(breaker.container(), filename + ".fastg", filename + ".fasta");
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
    } else if (io::SequencingLibraryBase::IsContigLib(type)) {
        return cfg::get().pe_params.long_reads.contigs.filtering;
    }
    return cfg::get().pe_params.long_reads.single_reads.filtering;
}

double GetSingleReadsWeightPriorityThreshold(const io::LibraryType& type) {
    if (type == io::LibraryType::PacBioReads || type == io::LibraryType::SangerReads || type == io::LibraryType::NanoporeReads) {
        return cfg::get().pe_params.long_reads.pacbio_reads.weight_priority;
    } else if (io::SequencingLibraryBase::IsContigLib(type)) {
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
    } else if (io::SequencingLibraryBase::IsContigLib(type)) {
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

inline void FinalizePaths(PathContainer& paths, GraphCoverageMap& cover_map, size_t max_overlap, bool mate_pairs = false) {
    ContigWriter writer(cover_map.graph());
    PathExtendResolver resolver(cover_map.graph());

    resolver.removeOverlaps(paths, cover_map, max_overlap, cfg::get().pe_params.param_set.remove_overlaps, cfg::get().pe_params.param_set.cut_all_overlaps);
    if (mate_pairs) {
        resolver.RemoveMatePairEnds(paths, max_overlap);
    }
    if (cfg::get().avoid_rc_connections) {
        paths.FilterInterstandBulges();
    }
    paths.FilterEmptyPaths();
    if (!mate_pairs) {
        resolver.addUncoveredEdges(paths, cover_map);
    }
    paths.SortByLength();
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
    MPStage,
    FinalizingPEStage,
    Scaffold2015,
};

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
    bool is_mp = lib.type() == io::LibraryType::MatePairs;
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

    ContigWriter writer(gp.g);
    const auto& lib = cfg::get().ds.reads[lib_index];
    shared_ptr<ExtensionChooser> longReadEC =
            make_shared<LongReadsExtensionChooser>(gp.g, paths, GetSingleReadsFilteringThreshold(lib.type()),
                                                   GetSingleReadsWeightPriorityThreshold(lib.type()),
                                                   GetSingleReadsUniqueEdgePriorityThreshold(lib.type()));

    size_t max_repeat_length = std::max(10000ul, lib.data().read_length);
    INFO("max_repeat_length set to " << max_repeat_length);
    return make_shared<SimpleExtender>(gp, cov_map, longReadEC, max_repeat_length,  
            pset.loop_removal.max_loops, true, UseCoverageResolverForSingleReads(lib.type()));
}

inline shared_ptr<SimpleExtender> MakeLongEdgePEExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                                         size_t lib_index, const pe_config::ParamSetT& pset, bool investigate_loops) {
    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.clustered_indices, lib_index);
    SetSingleThresholdForLib(lib, pset, cfg::get().ds.reads[lib_index].data().pi_threshold);

    shared_ptr<WeightCounter> wc = make_shared<PathCoverWeightCounter>(gp.g, lib, pset.normalize_weight);
    shared_ptr<ExtensionChooser> extension = make_shared<LongEdgeExtensionChooser>(gp.g, wc, GetWeightThreshold(lib, pset), GetPriorityCoeff(lib, pset));
    return make_shared<SimpleExtender>(gp, cov_map, extension, lib->GetISMax(), pset.loop_removal.max_loops, investigate_loops, false);
}

inline shared_ptr<SimpleExtender> MakeMetaExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                       size_t lib_index, const pe_config::ParamSetT& pset, bool investigate_loops) {
    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.clustered_indices, lib_index);

    shared_ptr<WeightCounter> wc = make_shared<MetagenomicWeightCounter>(gp.g, lib, /*read_length*/cfg::get().ds.RL(), 
                            /*normalized_threshold*/ 0.3, /*raw_threshold*/ 3, /*estimation_edge_length*/ 500);
    shared_ptr<SimpleExtensionChooser> extension = make_shared<SimpleExtensionChooser>(gp.g, wc, GetWeightThreshold(lib, pset), GetPriorityCoeff(lib, pset));
    return make_shared<SimpleExtender>(gp, cov_map, extension, lib->GetISMax(), pset.loop_removal.max_loops, investigate_loops, false);
}

inline shared_ptr<SimpleExtender> MakePEExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                       size_t lib_index, const pe_config::ParamSetT& pset, bool investigate_loops) {
    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.clustered_indices, lib_index);
    SetSingleThresholdForLib(lib, pset, cfg::get().ds.reads[lib_index].data().pi_threshold);

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
    auto gap_joiner = std::make_shared<HammingGapJoiner>(gp.g, pset.scaffolder_options.min_gap_score,
                                                 int(math::round((double) gp.g.k() - var_coeff * (double) lib->GetIsVar())),
                                                 (int) (pset.scaffolder_options.max_can_overlap * (double) gp.g.k()),
                                                 pset.scaffolder_options.short_overlap,
                                                 (int) 2 * cfg::get().ds.RL(), pset.scaffolder_options.artificial_gap,
                                                 cfg::get().pe_params.param_set.scaffolder_options.use_old_score);
    return make_shared<ScaffoldingPathExtender>(gp, cov_map, scaff_chooser, gap_joiner, lib->GetISMax(), pset.loop_removal.max_loops, false);
}


inline shared_ptr<PathExtender> MakeScaffolding2015Extender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                                        size_t lib_index, const pe_config::ParamSetT& pset, shared_ptr<ScaffoldingUniqueEdgeStorage> storage) {
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
    double prior_coef = GetPriorityCoeff(lib, pset);
    double var_coeff = 3.0;
    DEBUG("here creating extchooser");
    //TODO: 2 is relative weight cutoff, to config!
    auto scaff_chooser = std::make_shared<ExtensionChooser2015>(gp.g, counter, prior_coef, var_coeff, storage, 2);
    //);
    auto gap_joiner = std::make_shared<HammingGapJoiner>(gp.g, pset.scaffolder_options.min_gap_score,
                                                         int(math::round((double) gp.g.k() - var_coeff * (double) lib->GetIsVar())),
                                                         (int) (pset.scaffolder_options.max_can_overlap * (double) gp.g.k()),
                                                         pset.scaffolder_options.short_overlap,
                                                         (int) 2 * cfg::get().ds.RL(), pset.scaffolder_options.artificial_gap,
                                                         cfg::get().pe_params.param_set.scaffolder_options.use_old_score);
    return make_shared<ScaffoldingPathExtender>(gp, cov_map, scaff_chooser, gap_joiner, lib->GetISMax(), pset.loop_removal.max_loops, false , false);
}


inline shared_ptr<SimpleExtender> MakeMPExtender(const conj_graph_pack& gp, const GraphCoverageMap& cov_map, const PathContainer& paths,
                                       size_t lib_index, const pe_config::ParamSetT& pset) {

    shared_ptr<PairedInfoLibrary> lib = MakeNewLib(gp.g, gp.paired_indices, lib_index);
    size_t max_number_of_paths_to_search = GetNumberMPPaths(gp.g);
    DEBUG("max number of mp paths " << max_number_of_paths_to_search);

    shared_ptr<MatePairExtensionChooser> chooser = make_shared<MatePairExtensionChooser>(gp.g, lib, paths, max_number_of_paths_to_search);
    return make_shared<SimpleExtender>(gp, cov_map, chooser, lib->GetISMax(), pset.loop_removal.mp_max_loops, true, false);
}


inline bool InsertSizeCompare(const shared_ptr<PairedInfoLibrary> lib1,
                              const shared_ptr<PairedInfoLibrary> lib2) {
    return lib1->GetISMax() < lib2->GetISMax();
}

inline vector<shared_ptr<PathExtender> > MakeAllScaffoldingExtenders2015(PathExtendStage stage,
                                                                         const conj_graph_pack &gp,
                                                                         const GraphCoverageMap &cov_map,
                                                                         const pe_config::ParamSetT &pset) {
    VERIFY_MSG(false, "Do not use MakeAllScaffoldingExtenders2015");
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp, cfg::get().pe_params.scaffolding2015.min_unique_length, cfg::get().pe_params.scaffolding2015.unique_coverage_variation);
    auto storage = std::make_shared<ScaffoldingUniqueEdgeStorage>();

    unique_edge_analyzer.FillUniqueEdgeStorage(*storage);
//    GenomeConsistenceChecker genome_checker (gp, *storage, 500, 0.2);
//    genome_checker.SpellGenome();
    vector<shared_ptr<PathExtender> > result;
    DEBUG(cfg::get().ds.reads.lib_count());
    for (io::LibraryType lt : io::LibraryPriotity) {
        DEBUG("strt");
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            const auto &lib = cfg::get().ds.reads[i];                        
            if (lib.type() != lt)
                continue;
            if (IsForMPExtender(lib) || IsForPEExtender(lib)) {
                DEBUG("pushing");
                result.push_back(MakeScaffolding2015Extender(gp, cov_map, i, pset, storage));
            } else {
                DEBUG("sinlge read only")
            }
        }
    }
    return result;
}

template<typename Base, typename T>
inline bool instanceof(const T *ptr) {
    return dynamic_cast<const Base*>(ptr) != nullptr;
}

inline void PrintExtenders(vector<shared_ptr<PathExtender> >& extenders) {
    INFO("Extenders in vector:");
    for(size_t i = 0; i < extenders.size(); ++i) {
        string type = typeid(*extenders[i]).name();
        INFO("Extender #i" << type);
        if (instanceof<SimpleExtender>(extenders[i].get())) {
            auto ec = ((SimpleExtender *) extenders[i].get())->GetExtensionChooser();
            string chooser_type = typeid(*ec).name();
            INFO("    Extender #i" << chooser_type);
        }
        else if (instanceof<ScaffoldingPathExtender>(extenders[i].get())) {
            auto ec = ((ScaffoldingPathExtender *) extenders[i].get())->GetExtensionChooser();
            string chooser_type = typeid(*ec).name();
            INFO("    Extender #i" << chooser_type);
        }
    }
}

inline vector<shared_ptr<PathExtender> > MakeAllExtenders(PathExtendStage stage, const conj_graph_pack& gp, const GraphCoverageMap& cov_map,
                                            const pe_config::ParamSetT& pset, shared_ptr<ScaffoldingUniqueEdgeStorage> storage, const PathContainer& paths_for_mp = PathContainer()) {

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

            if (IsForSingleReadExtender(lib)  && pset.sm != sm_2015) {
                result.push_back(MakeLongReadsExtender(gp, cov_map, i, pset));
                ++single_read_libs;
            }
            if (IsForPEExtender(lib) && stage == PathExtendStage::PEStage) {
                ++pe_libs;
                if (pset.sm == sm_old_pe_2015 || pset.sm == sm_old || pset.sm == sm_combined) {
                    if (cfg::get().ds.moleculo)
                        pes.push_back(MakeLongEdgePEExtender(gp, cov_map, i, pset, false));
                    pes.push_back(MakePEExtender(gp, cov_map, i, pset, false));
                }
            }
            if (IsForShortLoopExtender(lib) && (pset.sm == sm_old_pe_2015 || pset.sm == sm_old || pset.sm == sm_combined)) {
                pe_loops.push_back(MakePEExtender(gp, cov_map, i, pset, true));
            }
            if (IsForPEExtender(lib) && pset.sm == sm_2015) {
                pes2015.push_back(MakeScaffolding2015Extender(gp, cov_map, i, pset, storage));
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
            if (IsForMPExtender(lib) && stage == PathExtendStage::MPStage) {
                ++mp_libs;
                if (pset.sm == sm_old || pset.sm == sm_combined) {
                    mps.push_back(MakeMPExtender(gp, cov_map, paths_for_mp, i, pset));
                }
                if (pset.sm == sm_old_pe_2015 || pset.sm == sm_2015 || pset.sm == sm_combined) {
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
    PrintExtenders(result);
    return result;
}

size_t FindOverlapLenForStage(PathExtendStage stage) {
    size_t res = 0;
    for (const auto& lib : cfg::get().ds.reads) {
        if (IsForPEExtender(lib) && stage == PathExtendStage::PEStage) {
            res = max(res, (size_t) lib.data().insert_size_right_quantile);
        } else if (IsForShortLoopExtender(lib)) {
            res = max(res, (size_t) lib.data().insert_size_right_quantile);
        } else if (IsForMPExtender(lib) && stage == PathExtendStage::MPStage) {
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

inline void FinalizeUniquenessPaths(){
//TODO::compilation placeholder
}


inline void ResolveRepeatsPe(conj_graph_pack& gp,
        const std::string& output_dir,
        const std::string& contigs_name,
        bool traversLoops,
        boost::optional<std::string> broken_contigs) {

    INFO("ExSPAnder repeat resolving tool started");


    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp, cfg::get().pe_params.scaffolding2015.min_unique_length, cfg::get().pe_params.scaffolding2015.unique_coverage_variation);
    auto storage = std::make_shared<ScaffoldingUniqueEdgeStorage>();
    unique_edge_analyzer.FillUniqueEdgeStorage(*storage);
    auto sc_mode = cfg::get().pe_params.param_set.sm;

    make_dir(output_dir);
    make_dir(GetEtcDir(output_dir));
    const pe_config::ParamSetT& pset = cfg::get().pe_params.param_set;

    ContigWriter writer(gp.g);

//make pe + long reads extenders
    GraphCoverageMap cover_map(gp.g);
    INFO("SUBSTAGE = paired-end libraries")
    PathExtendStage exspander_stage = PathExtendStage::PEStage;
    vector<shared_ptr<PathExtender> > all_libs = MakeAllExtenders(exspander_stage, gp, cover_map, pset, storage);


    size_t max_over = max(FindOverlapLenForStage(exspander_stage), gp.g.k() + 100);
    shared_ptr<CompositeExtender> mainPE = make_shared<CompositeExtender>(gp.g, cover_map, all_libs, max_over, storage);

//extend pe + long reads
    PathExtendResolver resolver(gp.g);
    auto seeds = resolver.makeSimpleSeeds();
    DebugOutputPaths(writer, gp, output_dir, seeds, "init_paths");
    seeds.SortByLength();
    INFO("Growing paths using paired-end and long single reads");
    auto paths = resolver.extendSeeds(seeds, *mainPE);
    paths.SortByLength();
    DebugOutputPaths(writer, gp, output_dir, paths, "pe_overlaped_paths");

    PathContainer clone_paths;
    GraphCoverageMap clone_map(gp.g);
    bool mp_exist = MPLibsExist();

    if (mp_exist) {
        ClonePathContainer(paths, clone_paths, clone_map);
    }

    if (!(sc_mode == sm_old_pe_2015 || sc_mode == sm_2015 || sc_mode == sm_combined))
        FinalizePaths(paths, cover_map, max_over);
    if (broken_contigs.is_initialized()) {
        OutputBrokenScaffolds(paths, (int) gp.g.k(), writer,
                              output_dir + (mp_exist ? "pe_contigs" : broken_contigs.get()));
    }
    writer.WritePathsToFASTG(paths, GetEtcDir(output_dir) + "pe_before_traversal.fastg", GetEtcDir(output_dir) + "pe_before_traversal.fasta");
    DebugOutputPaths(writer, gp, output_dir, paths, "before_traverse_pe");
    if (traversLoops) {
        TraverseLoops(paths, cover_map, mainPE);
    }
    DebugOutputPaths(writer, gp, output_dir, paths, (mp_exist ? "final_pe_paths" : "final_paths"));
    writer.WritePathsToFASTG(paths,
                             output_dir + (mp_exist ? "pe_scaffolds" : contigs_name) + ".fastg",
                             output_dir + (mp_exist ? "pe_scaffolds" : contigs_name) + ".fasta" );

    cover_map.Clear();
    paths.DeleteAllPaths();
    if (!mp_exist) {
        return;
    }

//MP
    DebugOutputPaths(writer, gp, output_dir, clone_paths, "before_mp_paths");

    INFO("SUBSTAGE = mate-pair libraries ")
    exspander_stage = PathExtendStage::MPStage;
    all_libs.clear();
    all_libs = MakeAllExtenders(exspander_stage, gp, clone_map, pset, storage, clone_paths);
    max_over = FindOverlapLenForStage(exspander_stage);
    shared_ptr<CompositeExtender> mp_main_pe = make_shared<CompositeExtender>(gp.g, clone_map, all_libs, max_over, storage);

    INFO("Growing paths using mate-pairs");
    auto mp_paths = resolver.extendSeeds(clone_paths, *mp_main_pe);
    if (!(sc_mode == sm_old_pe_2015 || sc_mode == sm_2015 || sc_mode == sm_combined))
        FinalizePaths(mp_paths, clone_map, max_over, true);
    DebugOutputPaths(writer, gp, output_dir, mp_paths, "mp_final_paths");
    writer.WritePathsToFASTG(mp_paths, GetEtcDir(output_dir) + "mp_prefinal.fastg", GetEtcDir(output_dir) + "mp_prefinal.fasta");

    DEBUG("Paths are grown with mate-pairs");
    if (cfg::get().pe_params.debug_output) {
        writer.writePaths(mp_paths, output_dir + "mp_paths.fasta");
    }
//MP end

//pe again
    INFO("SUBSTAGE = polishing paths")
    exspander_stage = PathExtendStage::FinalizingPEStage;
    all_libs.clear();
    all_libs = MakeAllExtenders(exspander_stage, gp, cover_map, pset, storage);
    max_over = FindOverlapLenForStage(exspander_stage);
    shared_ptr<CompositeExtender> last_extender = make_shared<CompositeExtender>(gp.g, clone_map, all_libs, max_over, storage);

    auto last_paths = resolver.extendSeeds(mp_paths, *last_extender);
    if (!(sc_mode == sm_old_pe_2015 || sc_mode == sm_2015 || sc_mode == sm_combined))
        FinalizePaths(last_paths, clone_map, max_over);
    writer.WritePathsToFASTG(last_paths, GetEtcDir(output_dir) + "mp_before_traversal.fastg", GetEtcDir(output_dir) + "mp_before_traversal.fasta");
    DebugOutputPaths(writer, gp, output_dir, last_paths, "before_traverse_mp");
    TraverseLoops(last_paths, clone_map, last_extender);

//result
    if (broken_contigs.is_initialized()) {
        OutputBrokenScaffolds(last_paths, (int) gp.g.k(), writer, output_dir + broken_contigs.get());
    }
    DebugOutputPaths(writer, gp, output_dir, last_paths, "last_paths");
    writer.WritePathsToFASTG(last_paths, output_dir + contigs_name + ".fastg", output_dir + contigs_name + ".fasta" , gp);

    //FinalizeUniquenessPaths();

    last_paths.DeleteAllPaths();
    seeds.DeleteAllPaths();
    mp_paths.DeleteAllPaths();
    clone_paths.DeleteAllPaths();

    INFO("ExSPAnder repeat resolving tool finished");
}
/*
inline void ScaffoldAll2015(conj_graph_pack& gp,
                             const std::string& output_dir,
                             const std::string& contigs_name,
                             bool traversLoops,
                             boost::optional<std::string> broken_contigs,
                             bool use_auto_threshold = true) {

    INFO("Scaffolding 2015 started");

    make_dir(output_dir);
    make_dir(GetEtcDir(output_dir));
    const pe_config::ParamSetT& pset = cfg::get().pe_params.param_set;

    ContigWriter writer(gp.g);

//make pe + long reads extenders
    GraphCoverageMap cover_map(gp.g);
    INFO("SUBSTAGE = paired-end libraries")
    PathExtendStage exspander_stage = PathExtendStage::Scaffold2015;
    vector<shared_ptr<PathExtender> > all_libs = MakeAllScaffoldingExtenders2015(exspander_stage, gp, cover_map, pset);
    size_t max_over = max(FindOverlapLenForStage(exspander_stage), gp.g.k() + 100);
    shared_ptr<CompositeExtender> main_extender = make_shared<CompositeExtender>(gp.g, cover_map, all_libs, max_over);

//extend pe + long reads
    PathExtendResolver resolver(gp.g);
    auto seeds = resolver.makeSimpleSeeds();
    DebugOutputPaths(writer, gp, output_dir, seeds, "init_paths");
    seeds.SortByLength();
    INFO("Growing paths using paired-end and long single reads");
    auto paths = resolver.extendSeeds(seeds, *main_extender);
    paths.SortByLength();
    DebugOutputPaths(writer, gp, output_dir, paths, "pe_overlaped_paths");

    PathContainer clone_paths;
    GraphCoverageMap clone_map(gp.g);
    bool mp_exist = MPLibsExist();

    FinalizeUniquenessPaths();

    writer.WritePathsToFASTG(paths, GetEtcDir(output_dir) + "scaf_before_traversal.fastg", GetEtcDir(output_dir) + "scaf_before_traversal.fasta");
    DebugOutputPaths(writer, gp, output_dir, paths, "before_traverse_scaf");
    if (traversLoops) {
        TraverseLoops(paths, cover_map, main_extender);
    }
    DebugOutputPaths(writer, gp, output_dir, paths, (mp_exist ? "final_pe_paths" : "final_paths"));
    writer.WritePathsToFASTG(paths,
                             output_dir + ("scaffolds2015") + ".fastg",
                             output_dir + ("scaffolds2015") + ".fasta" , gp);

    cover_map.Clear();

    paths.DeleteAllPaths();
    seeds.DeleteAllPaths();
    clone_paths.DeleteAllPaths();

    INFO("ExSPAnder repeat resolving tool finished");
}

*/
} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
