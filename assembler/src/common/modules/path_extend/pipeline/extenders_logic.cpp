//
// Created by andrey on 14.11.16.
//

#include "extenders_logic.hpp"


namespace path_extend {

using namespace debruijn_graph;

void ExtendersGenerator::AddPathsToContainer(const std::vector<PathInfo<Graph>> &paths,
                                             size_t size_threshold,
                                             PathContainer &result) const {
    for (const auto &path: paths) {
        const auto &edges = path.getPath();
        if (edges.size() <= size_threshold)
            continue;

        BidirectionalPath *new_path = new BidirectionalPath(gp_.g, edges);
        BidirectionalPath *conj_path = new BidirectionalPath(new_path->Conjugate());
        new_path->SetWeight((float) path.getWeight());
        conj_path->SetWeight((float) path.getWeight());
        result.AddPair(new_path, conj_path);
    }
    DEBUG("Long reads paths " << result.size() << " == ");
}

shared_ptr<ExtensionChooser> ExtendersGenerator::MakeLongReadsExtensionChooser(const config::dataset::Library &lib,
                                                                               size_t lib_index) const {
    PathContainer paths;
    AddPathsToContainer(gp_.single_long_reads[lib_index].GetAllPaths(), 1, paths);

    auto long_reads_config = support_.GetLongReadsConfig(lib.type());
    return make_shared<LongReadsExtensionChooser>(gp_.g, paths, long_reads_config.filtering,
                                                  long_reads_config.weight_priority,
                                                  long_reads_config.unique_edge_priority,
                                                  long_reads_config.min_significant_overlap,
                                                  params_.pset.extension_options.max_repeat_length,
                                                  params_.uneven_depth);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeLongReadsExtender(size_t lib_index) const {
    const auto &lib = dataset_info_.reads[lib_index];
    //TODO params
    //FIXME what would you suggest?
    size_t resolvable_repeat_length_bound = 10000ul;
    if (!lib.is_contig_lib()) {
        resolvable_repeat_length_bound = std::max(resolvable_repeat_length_bound, lib.data().read_length);
    }
    INFO("resolvable_repeat_length_bound set to " << resolvable_repeat_length_bound);


    auto long_read_ec = MakeLongReadsExtensionChooser(lib, lib_index);
    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       long_read_ec,
                                       resolvable_repeat_length_bound,
                                       params_.pset.loop_removal.max_loops,
                                       true, /* investigate short loops */
                                       support_.UseCoverageResolverForSingleReads(lib.type()));
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeLongEdgePEExtender(size_t lib_index,
                                                                      bool investigate_loops) const {
    const auto &lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
    VERIFY_MSG(math::ge(lib.data().pi_threshold, 0.0), "PI threshold should be set");
    support_.SetSingleThresholdForLib(paired_lib, params_.pset, lib.data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<WeightCounter> wc =
        make_shared<PathCoverWeightCounter>(gp_.g, paired_lib, params_.pset.normalize_weight);
    auto opts = support_.GetExtensionOpts(paired_lib, params_.pset);
    shared_ptr<ExtensionChooser> extension =
        make_shared<LongEdgeExtensionChooser>(gp_.g, wc,
                                              opts.weight_threshold,
                                              opts.priority_coeff);

    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       extension,
                                       paired_lib->GetISMax(),
                                       params_.pset.loop_removal.max_loops,
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/);
}

shared_ptr<SimpleExtensionChooser> ExtendersGenerator::MakeMetaExtensionChooser(shared_ptr<PairedInfoLibrary> lib,
                                                                                size_t read_length) const {
    VERIFY(params_.mode == config::pipeline_type::meta);
    VERIFY(!lib->IsMp());
    shared_ptr<WeightCounter> wc = make_shared<MetagenomicWeightCounter>(gp_.g,
                                                                         lib,
                                                                         read_length, //read_length
                                                                         0.3, //normalized_threshold
                                                                         3, //raw_threshold
                                                                         0 /*estimation_edge_length*/ );
    return make_shared<SimpleExtensionChooser>(gp_.g, wc,
                                               params_.pset.extension_options.weight_threshold,
                                               params_.pset.extension_options.priority_coeff);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeMetaExtender(size_t lib_index,
                                                                bool investigate_loops) const {

    const auto &lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);

    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       MakeMetaExtensionChooser(paired_lib, dataset_info_.RL()),
                                       paired_lib->GetISMax(),
                                       params_.pset.loop_removal.max_loops,
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/);
}


shared_ptr<GapJoiner> ExtendersGenerator::MakeGapJoiners(shared_ptr<PairedInfoLibrary> paired_lib) const {
    const auto &pset = params_.pset;

    vector<shared_ptr<GapJoiner>> joiners;
    if (params_.pset.scaffolder_options.use_la_gap_joiner)
        joiners.push_back(std::make_shared<LAGapJoiner>(gp_.g, pset.scaffolder_options.min_overlap_length,
                                                        pset.scaffolder_options.flank_multiplication_coefficient,
                                                        pset.scaffolder_options.flank_addition_coefficient));


    joiners.push_back(std::make_shared<HammingGapJoiner>(gp_.g,
                                                         pset.scaffolder_options.min_gap_score,
                                                         pset.scaffolder_options.short_overlap,
                                                         (int) pset.scaffolder_options.basic_overlap_coeff
                                                             * dataset_info_.RL()));

    return std::make_shared<CompositeGapJoiner>(gp_.g,
                                                joiners,
                                                size_t(pset.scaffolder_options.max_can_overlap
                                                           * (double) gp_.g.k()), /* may overlap threshold */
                                                int(math::round((double) gp_.g.k()
                                                                    - pset.scaffolder_options.var_coeff
                                                                        * (double) paired_lib->GetIsVar())),  /* must overlap threshold */
                                                pset.scaffolder_options.artificial_gap);

}

shared_ptr<PathExtender> ExtendersGenerator::MakeScaffoldingExtender(size_t lib_index) const {

    const auto &lib = dataset_info_.reads[lib_index];
    const auto &pset = params_.pset;
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.scaffolding_indices[lib_index]);

    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(gp_.g, paired_lib);

    auto scaff_chooser = std::make_shared<ScaffoldingExtensionChooser>(gp_.g, counter,
                                                                       pset.scaffolder_options.cl_threshold,
                                                                       pset.scaffolder_options.var_coeff);

    return make_shared<ScaffoldingPathExtender>(gp_, cover_map_, scaff_chooser,
                                                MakeGapJoiners(paired_lib),
                                                paired_lib->GetISMax(),
                                                pset.loop_removal.max_loops,
                                                false, /* investigate short loops */
                                                params_.avoid_rc_connections);
}

shared_ptr<PathExtender> ExtendersGenerator::MakeRNAScaffoldingExtender(size_t lib_index) const {

    const auto &lib = dataset_info_.reads[lib_index];
    const auto &pset = params_.pset;
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.paired_indices[lib_index]);

    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(gp_.g, paired_lib);

    auto scaff_chooser = std::make_shared<ScaffoldingExtensionChooser>(gp_.g,
                                                                       counter,
                                                                       pset.scaffolder_options.cutoff,
                                                                       pset.scaffolder_options.var_coeff);
    auto scaff_chooser2 = std::make_shared<ScaffoldingExtensionChooser>(gp_.g,
                                                                        counter,
                                                                        pset.scaffolder_options.hard_cutoff,
                                                                        pset.scaffolder_options.var_coeff);


    VERIFY(pset.scaffolder_options.min_overlap_for_rna_scaffolding.is_initialized());
    return make_shared<RNAScaffoldingPathExtender>(gp_, cover_map_,
                                                   scaff_chooser,
                                                   scaff_chooser2,
                                                   MakeGapJoiners(paired_lib),
                                                   paired_lib->GetISMax(),
                                                   pset.loop_removal.max_loops,
                                                   false  /* investigate short loops */,
                                                   *pset.scaffolder_options.min_overlap_for_rna_scaffolding);
}

//FIXME rename!
shared_ptr<PathExtender> ExtendersGenerator::MakeScaffolding2015Extender(
    size_t lib_index,
    const ScaffoldingUniqueEdgeStorage &storage) const {

    const auto &lib = dataset_info_.reads[lib_index];
    const auto &pset = params_.pset;
    shared_ptr<PairedInfoLibrary> paired_lib;
    INFO("Creating Scaffolding 2015 extender for lib #" << lib_index);

    //TODO:: temporary solution
    //weird check!
    //ANSWER: DimaA heritage, will fix
    if (gp_.paired_indices[lib_index].size() > gp_.clustered_indices[lib_index].size()) {
        INFO("Paired unclustered indices not empty, using them");
        paired_lib = MakeNewLib(gp_.g, lib, gp_.paired_indices[lib_index]);
    } else if (gp_.clustered_indices[lib_index].size() != 0) {
        INFO("clustered indices not empty, using them");
        paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
    } else {
        ERROR("All paired indices are empty!");
    }

    //TODO::was copypasted from MakeScaffoldingExtender, refactor 2015 extension chooser
    DEBUG("creating extchooser");
    shared_ptr<ConnectionCondition>
        condition = make_shared<PairedLibConnectionCondition>(gp_.g, paired_lib, lib_index, 0);
    auto scaff_chooser = std::make_shared<ExtensionChooser2015>(gp_.g,
                                                                nullptr,
                                                                condition,
                                                                storage,
                                                                pset.scaffolder_options.cl_threshold,
                                                                pset.scaffolder_options.var_coeff,
                                                                pset.scaffolding2015.relative_weight_cutoff,
                                                                gp_.g.size()
                                                                    <= params_.pset.scaffolding2015.graph_connectivity_max_edges);

    //TODO: should we switch to normal gap joiner?
    auto gap_joiner = std::make_shared<HammingGapJoiner>(gp_.g, pset.scaffolder_options.min_gap_score,
                                                         pset.scaffolder_options.short_overlap,
                                                         (int) pset.scaffolder_options.basic_overlap_coeff
                                                             * dataset_info_.RL());

    return make_shared<ScaffoldingPathExtender>(gp_, cover_map_,
                                                scaff_chooser,
                                                gap_joiner,
                                                paired_lib->GetISMax(),
                                                pset.loop_removal.max_loops,
                                                false, /* investigate short loops */
                                                params_.avoid_rc_connections,
                                                false /* jump only from tips */);
}

//FIXME to DimaM: magic constants and -1ul need comments
shared_ptr<SimpleExtender> ExtendersGenerator::MakeCoordCoverageExtender(size_t lib_index) const {

    const auto &lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);

    CoverageAwareIdealInfoProvider provider(gp_.g, paired_lib, -1ul, 0);
    auto coord_chooser = make_shared<CoordinatedCoverageExtensionChooser>(gp_.g, provider,
                                                                          params_.pset.coordinated_coverage.max_edge_length_in_repeat,
                                                                          params_.pset.coordinated_coverage.delta,
                                                                          params_.pset.coordinated_coverage.min_path_len);
    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(gp_.g, paired_lib, false);
    auto chooser = make_shared<JointExtensionChooser>(gp_.g,
                                                      make_shared<TrivialExtensionChooserWithPI>(gp_.g,
                                                                                                 counter,
                                                                                                 1.5),
                                                      coord_chooser);
    return make_shared<SimpleExtender>(gp_,
                                       cover_map_,
                                       chooser,
                                       -1ul,
                                       params_.pset.loop_removal.mp_max_loops,
                                       false,
                                       false);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeRNAExtender(size_t lib_index, bool investigate_loops) const {

    const auto &lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);

    support_.SetSingleThresholdForLib(paired_lib, params_.pset, lib.data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<WeightCounter>
        wc = make_shared<PathCoverWeightCounter>(gp_.g, paired_lib, params_.pset.normalize_weight);
    auto opts = support_.GetExtensionOpts(paired_lib, params_.pset);
    shared_ptr<RNAExtensionChooser> extension =
        make_shared<RNAExtensionChooser>(gp_.g, wc,
                                         opts.weight_threshold,
                                         opts.priority_coeff);

    return make_shared<MultiExtender>(gp_, cover_map_,
                                      extension,
                                      paired_lib->GetISMax(),
                                      params_.pset.loop_removal.max_loops,
                                      investigate_loops,
                                      false /*use short loop coverage resolver*/);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakePEExtender(
    size_t lib_index,
    bool investigate_loops) const {

    const auto &lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
    VERIFY_MSG(math::ge(lib.data().pi_threshold, 0.0), "PI threshold should be set");
    support_.SetSingleThresholdForLib(paired_lib, params_.pset, lib.data().pi_threshold);
    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<WeightCounter>
        wc = make_shared<PathCoverWeightCounter>(gp_.g, paired_lib, params_.pset.normalize_weight);
    auto opts = support_.GetExtensionOpts(paired_lib, params_.pset);
    auto extension = make_shared<SimpleExtensionChooser>(gp_.g, wc,
                                                         opts.weight_threshold,
                                                         opts.priority_coeff);

    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       extension,
                                       paired_lib->GetISMax(),
                                       params_.pset.loop_removal.max_loops,
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/);
}


void ExtendersGenerator::PrintExtenders(Extenders &extenders) const {
    DEBUG("Extenders in vector:");
    for (size_t i = 0; i < extenders.size(); ++i) {
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

ExtendersGenerator::Extenders ExtendersGenerator::MakeMPExtenders(const ScaffoldingUniqueEdgeStorage &storage) const {

    vector<shared_ptr<PathExtender> > result;
    size_t mp_libs = 0;

    for (io::LibraryType lt : io::LibraryPriotity) {
        for (size_t i = 0; i < dataset_info_.reads.lib_count(); ++i) {
            const auto &lib = dataset_info_.reads[i];
            if (lib.type() != lt)
                continue;

            if (lib.is_mate_pair()) {
                ++mp_libs;
                result.push_back(MakeScaffolding2015Extender(i, storage));
            }
        }
    }
    PrintExtenders(result);
    return result;
}

//FIXME why all parameters non const?
//To DimaA
ExtendersGenerator::Extenders ExtendersGenerator::MakePBScaffoldingExtenders(ScaffoldingUniqueEdgeStorage &unique_storage_pb,
                                                                             vector<PathContainer> &long_reads_paths,
                                                                             vector<shared_ptr<GraphCoverageMap>> &long_reads_cov_map) const {
    const auto &pset = params_.pset;
    //FIXME magic constants
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer_pb(gp_, 500, 0.5);
    Extenders result;
    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
        AddPathsToContainer(gp_.single_long_reads[lib_index].GetAllPaths(), 1, long_reads_paths[lib_index]);
        auto coverage_map = make_shared<GraphCoverageMap>(gp_.g, long_reads_paths[lib_index]);
        long_reads_cov_map.push_back(coverage_map);
    }
    INFO("Filling backbone edges for long reads scaffolding ...");
    if (params_.uneven_depth) {
        INFO("with long reads paths.");
        //TODO:: muiltiple libraries?
        for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
            if (dataset_info_.reads[lib_index].type() == io::LibraryType::TSLReads) {

                unique_edge_analyzer_pb.FillUniqueEdgesWithLongReads(long_reads_cov_map[lib_index],
                                                                     unique_storage_pb,
                                                                     support_.GetLongReadsConfig(dataset_info_.reads[lib_index].type()));
            }
        }
        INFO("removing fake unique wuth paired-end libs");
        for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
            if (dataset_info_.reads[lib_index].type() == io::LibraryType::PairedEnd) {
                unique_edge_analyzer_pb.ClearLongEdgesWithPairedLib(lib_index, unique_storage_pb);
            }
        }

    } else {
        INFO("with coverage.")
        unique_edge_analyzer_pb.FillUniqueEdgeStorage(unique_storage_pb);
    }
    INFO(unique_storage_pb.size() << " unique edges");

    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
        if (support_.IsForSingleReadExtender(dataset_info_.reads[lib_index])) {
            INFO("creating scaffolding extender for lib " << lib_index);
            shared_ptr<ConnectionCondition> condition = make_shared<LongReadsLibConnectionCondition>(gp_.g,
                                                                                                     lib_index, 2,
                                                                                                     *long_reads_cov_map[lib_index]);
            auto scaff_chooser = std::make_shared<ExtensionChooser2015>(gp_.g,
                                                                        nullptr,
                                                                        condition,
                                                                        unique_storage_pb,
                                                                        pset.scaffolder_options.cl_threshold,
                                                                        pset.scaffolder_options.var_coeff,
                                                                        pset.scaffolding2015.relative_weight_cutoff);

            auto gap_joiner = std::make_shared<HammingGapJoiner>(gp_.g, pset.scaffolder_options.min_gap_score,
                                                                 pset.scaffolder_options.short_overlap,
                                                                 (int) pset.scaffolder_options.basic_overlap_coeff *
                                                                     dataset_info_.RL());

            result.push_back(make_shared<ScaffoldingPathExtender>(gp_, cover_map_,
                                                                  scaff_chooser,
                                                                  gap_joiner,
                                                                  10000, /* insert size */
                                                                  pset.loop_removal.max_loops,
                                                                  false, /* investigate short loops */
                                                                  params_.avoid_rc_connections,
                                                                  false /* jump only from tips */));

        }
    }
    return result;
}

ExtendersGenerator::Extenders ExtendersGenerator::MakeBasicExtenders(const ScaffoldingUniqueEdgeStorage &storage) const {
    Extenders result;
    Extenders pes;
    Extenders pes2015;
    Extenders pe_loops;
    Extenders pe_scafs;

    size_t single_read_libs = 0;
    size_t pe_libs = 0;
    size_t scf_pe_libs = 0;

    const auto &pset = params_.pset;

    for (io::LibraryType lt : io::LibraryPriotity) {
        for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); ++lib_index) {
            const auto &lib = dataset_info_.reads[lib_index];
            if (lib.type() != lt)
                continue;

            //TODO: scaff2015 does not need any single read libs?
            if (support_.IsForSingleReadExtender(lib)) {
                result.push_back(MakeLongReadsExtender(lib_index));
                ++single_read_libs;
            }
            if (support_.IsForPEExtender(lib)) {
                ++pe_libs;
                if (IsOldPEEnabled(pset.sm)) {
                    if (params_.mode == config::pipeline_type::meta)
                        //TODO proper configuration via config
                        pes.push_back(MakeMetaExtender(lib_index, false));
                    else if (params_.mode == config::pipeline_type::moleculo)
                        pes.push_back(MakeLongEdgePEExtender(lib_index, false));
                    else if (pset.multi_path_extend) {
                        pes.push_back(MakePEExtender(lib_index, false));
                        //TODO: Exclude from polishing stage
                        pes.push_back(MakeRNAExtender(lib_index, false));
                    }
                    else
                        pes.push_back(MakePEExtender(lib_index, false));
                }
                else if (pset.sm == sm_2015) {
                    pes2015.push_back(MakeScaffolding2015Extender(lib_index, storage));
                }
            }
            //FIXME logic is very cryptic!
            if (support_.IsForShortLoopExtender(lib) && IsOldPEEnabled(pset.sm)) {
                if (params_.mode == config::pipeline_type::meta)
                    pes.push_back(MakeMetaExtender(lib_index, true));
                else
                    pe_loops.push_back(MakePEExtender(lib_index, true));

            }
            if (support_.IsForScaffoldingExtender(lib) && params_.use_scaffolder
                && pset.scaffolder_options.enabled) {
                ++scf_pe_libs;
                if (params_.mode == config::pipeline_type::rna) {
                    pe_scafs.push_back(MakeRNAScaffoldingExtender(lib_index));
                }
                else {
                    switch (pset.sm) {
                        case sm_old: {
                            pe_scafs.push_back(MakeScaffoldingExtender(lib_index));
                            break;
                        }
                        case sm_old_pe_2015: {
                            pe_scafs.push_back(MakeScaffoldingExtender(lib_index));
                            break;
                        }
                        case sm_combined: {
                            pe_scafs.push_back(MakeScaffoldingExtender(lib_index));
                            pe_scafs.push_back(MakeScaffolding2015Extender(lib_index, storage));
                            break;
                        }
                        default:
                            break;
                    }
                }
            }

        }

        //I added push_back_all method at some point, but not insist on using it
        result.insert(result.end(), pes.begin(), pes.end());
        result.insert(result.end(), pes2015.begin(), pes2015.end());
        result.insert(result.end(), pe_scafs.begin(), pe_scafs.end());
        //FIXME make scope variables
        pes.clear();
        pe_scafs.clear();
        pes2015.clear();
        //What about building a vector of triplets: <priority/lib#/extender> and sorting it?
        //ANSWER: good idea
    }

    result.insert(result.end(), pe_loops.begin(), pe_loops.end());
    pe_loops.clear();

    INFO("Using " << pe_libs << " paired-end " << support_.LibStr(pe_libs));
    INFO("Using " << scf_pe_libs << " paired-end scaffolding " << support_.LibStr(scf_pe_libs));
    INFO("Using " << single_read_libs << " single read " << support_.LibStr(single_read_libs));
    INFO("Scaffolder is " << (pset.scaffolder_options.enabled ? "on" : "off"));


    if (pset.use_coordinated_coverage) {
        INFO("Using additional coordinated coverage extender");
        result.push_back(MakeCoordCoverageExtender(0 /* lib index */));
    }

    PrintExtenders(result);
    return result;
}

}
