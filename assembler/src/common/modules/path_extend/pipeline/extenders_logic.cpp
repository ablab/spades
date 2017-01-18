//
// Created by andrey on 14.11.16.
//

#include "extenders_logic.hpp"
#include "modules/path_extend/scaffolder2015/extension_chooser2015.hpp"


namespace path_extend {

using namespace debruijn_graph;

shared_ptr<ExtensionChooser> ExtendersGenerator::MakeLongReadsExtensionChooser(size_t lib_index,
                                                                               const GraphCoverageMap &read_paths_cov_map) const {
    auto long_reads_config = support_.GetLongReadsConfig(dataset_info_.reads[lib_index].type());
    return make_shared<LongReadsExtensionChooser>(gp_.g, read_paths_cov_map,
                                                  long_reads_config.filtering,
                                                  long_reads_config.weight_priority,
                                                  long_reads_config.unique_edge_priority,
                                                  long_reads_config.min_significant_overlap,
                                                  params_.pset.extension_options.max_repeat_length,
                                                  params_.uneven_depth);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeLongReadsExtender(size_t lib_index,
                                                                     const GraphCoverageMap &read_paths_cov_map) const {
    const auto &lib = dataset_info_.reads[lib_index];
    //TODO params
    size_t resolvable_repeat_length_bound = 10000ul;
    if (!dataset_info_.reads[lib_index].is_contig_lib()) {
        resolvable_repeat_length_bound = std::max(resolvable_repeat_length_bound, lib.data().read_length);
    }
    INFO("resolvable_repeat_length_bound set to " << resolvable_repeat_length_bound);


    auto long_read_ec = MakeLongReadsExtensionChooser(lib_index, read_paths_cov_map);
    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       long_read_ec,
                                       resolvable_repeat_length_bound,
                                       true, /* investigate short loops */
                                       support_.UseCoverageResolverForSingleReads(lib.type()));
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeLongEdgePEExtender(size_t lib_index,
                                                                      bool investigate_loops) const {
    const auto &lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
    //INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<WeightCounter> wc =
        make_shared<PathCoverWeightCounter>(gp_.g, paired_lib,
                                            params_.pset.normalize_weight,
                                            support_.SingleThresholdForLib(params_.pset, lib.data().pi_threshold));
    auto opts = support_.GetExtensionOpts(paired_lib, params_.pset);
    shared_ptr<ExtensionChooser> extension =
        make_shared<LongEdgeExtensionChooser>(gp_.g, wc,
                                              opts.weight_threshold,
                                              opts.priority_coeff);

    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       extension,
                                       paired_lib->GetISMax(),
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/);
}

shared_ptr<GapJoiner> ExtendersGenerator::MakeGapJoiners(double is_variation) const {
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
                                                int(math::round(double(gp_.g.k())
                                                                    - pset.scaffolder_options.var_coeff
                                                                        * is_variation)),  /* must overlap threshold */
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
                                                MakeGapJoiners(paired_lib->GetIsVar()),
                                                paired_lib->GetISMax(),
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
                                                   MakeGapJoiners(paired_lib->GetIsVar()),
                                                   paired_lib->GetISMax(),
                                                   false  /* investigate short loops */,
                                                   *pset.scaffolder_options.min_overlap_for_rna_scaffolding);
}

shared_ptr<PathExtender> ExtendersGenerator::MakeMatePairScaffoldingExtender(
    size_t lib_index,
    const ScaffoldingUniqueEdgeStorage &storage) const {

    const auto &lib = dataset_info_.reads[lib_index];
    const auto &pset = params_.pset;
    shared_ptr<PairedInfoLibrary> paired_lib;
    INFO("Creating Scaffolding 2015 extender for lib #" << lib_index);

    //FIXME: DimaA
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

    return make_shared<ScaffoldingPathExtender>(gp_, cover_map_,
                                                scaff_chooser,
                                                MakeGapJoiners(paired_lib->GetIsVar()),
                                                paired_lib->GetISMax(),
                                                false, /* investigate short loops */
                                                params_.avoid_rc_connections,
                                                false /* jump only from tips */);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeCoordCoverageExtender(size_t lib_index) const {
    const auto& lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);

    auto provider = make_shared<CoverageAwareIdealInfoProvider>(gp_.g, paired_lib, dataset_info_.RL());

    auto meta_wc = make_shared<PathCoverWeightCounter>(gp_.g, paired_lib,
                                                       params_.pset.normalize_weight,
                                                       support_.SingleThresholdForLib(params_.pset, lib.data().pi_threshold),
                                                       provider);

    auto permissive_pi_chooser = make_shared<IdealBasedExtensionChooser>(gp_.g,
                                                                         meta_wc,
                                                                         params_.pset.extension_options.weight_threshold,
                                                                         params_.pset.extension_options.priority_coeff);

    auto coord_cov_chooser = make_shared<CoordinatedCoverageExtensionChooser>(gp_.g, *provider,
                                                                              params_.pset.coordinated_coverage.max_edge_length_in_repeat,
                                                                              params_.pset.coordinated_coverage.delta,
                                                                              params_.pset.coordinated_coverage.min_path_len);

    auto chooser = make_shared<JointExtensionChooser>(gp_.g, permissive_pi_chooser, coord_cov_chooser);

    return make_shared<SimpleExtender>(gp_, cover_map_, chooser,
                                       -1ul /* insert size is needed only for loop detection, which is not needed in this case */,
                                       false, /* investigate short loops */
                                       false /*use short loop coverage resolver*/);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeRNAExtender(size_t lib_index, bool investigate_loops) const {

    const auto &lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
//    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    auto cip = make_shared<CoverageAwareIdealInfoProvider>(gp_.g, paired_lib, dataset_info_.RL());
    shared_ptr<WeightCounter> wc =
        make_shared<PathCoverWeightCounter>(gp_.g, paired_lib, params_.pset.normalize_weight,
                                            support_.SingleThresholdForLib(params_.pset, lib.data().pi_threshold),
                                            cip);

    auto opts = support_.GetExtensionOpts(paired_lib, params_.pset);
    shared_ptr<RNAExtensionChooser> extension =
        make_shared<RNAExtensionChooser>(gp_.g, wc,
                                         opts.weight_threshold,
                                         opts.priority_coeff);

    return make_shared<MultiExtender>(gp_, cover_map_,
                                      extension,
                                      paired_lib->GetISMax(),
                                      investigate_loops,
                                      false /*use short loop coverage resolver*/);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakePEExtender(size_t lib_index, bool investigate_loops) const {
    const auto &lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
    VERIFY_MSG(!paired_lib->IsMp(), "Tried to create PE extender for MP library");
    auto opts = params_.pset.extension_options;
//    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<CoverageAwareIdealInfoProvider> iip = nullptr;
    if (opts.use_default_single_threshold) {
        if (params_.uneven_depth) {
            iip = make_shared<CoverageAwareIdealInfoProvider>(gp_.g, paired_lib, dataset_info_.RL());
        } else {
            double lib_cov = support_.EstimateLibCoverage(lib_index);
            INFO("Estimated coverage of library #" << lib_index << " is " << lib_cov);
            iip = make_shared<GlobalCoverageAwareIdealInfoProvider>(gp_.g, paired_lib, dataset_info_.RL(), lib_cov);
        }
    }
    auto wc = make_shared<PathCoverWeightCounter>(gp_.g, paired_lib, params_.pset.normalize_weight,
                                                  support_.SingleThresholdForLib(params_.pset, lib.data().pi_threshold),
                                                  iip);

    auto extension_chooser = make_shared<SimpleExtensionChooser>(gp_.g, wc,
                                                         opts.weight_threshold,
                                                         opts.priority_coeff);

    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       extension_chooser,
                                       paired_lib->GetISMax(),
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/);
}


void ExtendersGenerator::PrintExtenders(const Extenders &extenders) const {
    DEBUG("Extenders in vector:");
    for (const auto& extender : extenders) {
        //TODO: use polymorphism instead of RTTI
        auto ext_ptr = extender.get();
        DEBUG("Extender #i" << typeid(*ext_ptr).name());
        if (instanceof<SimpleExtender>(ext_ptr)) {
            auto ec = ((SimpleExtender *) ext_ptr)->GetExtensionChooser();
            auto ec_ptr = ec.get();
            DEBUG("    Extender #i" << typeid(*ec_ptr).name());
        }
        else if (instanceof<ScaffoldingPathExtender>(ext_ptr)) {
            auto ec = ((ScaffoldingPathExtender *) ext_ptr)->GetExtensionChooser();
            auto ec_ptr = ec.get();
            DEBUG("    Extender #i" << typeid(*ec_ptr).name());
        }
    }
}

Extenders ExtendersGenerator::MakeMPExtenders(const ScaffoldingUniqueEdgeStorage &storage) const {
    ExtenderTriplets result;

    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info_.reads[lib_index];

        if (lib.is_mate_pair()) {
            result.emplace_back(lib.type(), lib_index, MakeMatePairScaffoldingExtender(lib_index, storage));
        }
    }
    std::stable_sort(result.begin(), result.end());

    return ExtractExtenders(result);
}

Extenders ExtendersGenerator::MakePBScaffoldingExtenders(const ScaffoldingUniqueEdgeStorage &unique_storage_pb,
                                                         const vector<shared_ptr<GraphCoverageMap>> &long_reads_cov_map) const {
    const auto &pset = params_.pset;
    ExtenderTriplets result;

    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
        if (support_.IsForSingleReadScaffolder(dataset_info_.reads[lib_index])) {
            INFO("Creating scaffolding extender for lib " << lib_index);
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

            result.emplace_back(dataset_info_.reads[lib_index].type(),
                                lib_index,
                                make_shared<ScaffoldingPathExtender>(gp_, cover_map_,
                                                                     scaff_chooser,
                                                                     MakeGapJoiners(1000), /* "IS vatiation" */
                                                                     10000, /* insert size */
                                                                     false, /* investigate short loops */
                                                                     params_.avoid_rc_connections,
                                                                     false /* jump only from tips */));

        }
    }
    INFO("Using " << result.size() << " long reads scaffolding " << support_.LibStr(result.size()));
    std::stable_sort(result.begin(), result.end());

    return ExtractExtenders(result);
}


Extenders ExtendersGenerator::MakeCoverageExtenders() const {
    Extenders result;

    INFO("Using additional coordinated coverage extender");
    result.push_back(MakeCoordCoverageExtender(0 /* lib index */));

    return result;
}

Extenders ExtendersGenerator::MakeBasicExtenders(const ScaffoldingUniqueEdgeStorage &storage,
                                                 const vector<shared_ptr<GraphCoverageMap>> &long_reads_cov_map) const {
    ExtenderTriplets basic_extenders;
    ExtenderTriplets loop_resolving_extenders;
    ExtenderTriplets scaffolding_extenders;

    size_t single_read_libs = 0;
    size_t pe_libs = 0;
    size_t scf_pe_libs = 0;

    const auto &pset = params_.pset;

    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info_.reads[lib_index];

        //TODO: scaff2015 does not need any single read libs?
        if (support_.IsForSingleReadExtender(lib)) {
            basic_extenders.emplace_back(lib.type(), lib_index, MakeLongReadsExtender(lib_index, *long_reads_cov_map[lib_index]));
            ++single_read_libs;
        }
        if (support_.IsForPEExtender(lib)) {
            ++pe_libs;
            if (IsOldPEEnabled(pset.sm)) {
                if (params_.mode == config::pipeline_type::moleculo) {
                    basic_extenders.emplace_back(lib.type(), lib_index, MakeLongEdgePEExtender(lib_index, false));
                } else if (pset.multi_path_extend) {
                    basic_extenders.emplace_back(lib.type(), lib_index, MakePEExtender(lib_index, false));
                    basic_extenders.emplace_back(lib.type(), lib_index, MakeRNAExtender(lib_index, false));
                } else {
                    basic_extenders.emplace_back(lib.type(), lib_index, MakePEExtender(lib_index, false));
                }
            } else if (pset.sm == sm_2015) {
                basic_extenders.emplace_back(lib.type(), lib_index, MakeMatePairScaffoldingExtender(lib_index, storage));
            }
        }
        //TODO logic is very cryptic!
        if (support_.IsForShortLoopExtender(lib) && IsOldPEEnabled(pset.sm)) {
            loop_resolving_extenders.emplace_back(lib.type(), lib_index, MakePEExtender(lib_index, true));
            //TODO what about moleculo and rna here?
        }
        if (support_.IsForScaffoldingExtender(lib) && params_.use_scaffolder
            && pset.scaffolder_options.enabled) {
            ++scf_pe_libs;
            if (params_.mode == config::pipeline_type::rna) {
                scaffolding_extenders.emplace_back(lib.type(), lib_index, MakeRNAScaffoldingExtender(lib_index));
            } else {
                scaffolding_extenders.emplace_back(lib.type(), lib_index, MakeScaffoldingExtender(lib_index));
                if (pset.sm == sm_combined) {
                    scaffolding_extenders.emplace_back(lib.type(), lib_index, MakeMatePairScaffoldingExtender(lib_index, storage));
                }
            }
        }
    }

    std::stable_sort(basic_extenders.begin(), basic_extenders.end());
    std::stable_sort(scaffolding_extenders.begin(), scaffolding_extenders.end());
    std::stable_sort(loop_resolving_extenders.begin(), loop_resolving_extenders.end());

    Extenders result;
    push_back_all(result, ExtractExtenders(basic_extenders));
    push_back_all(result, ExtractExtenders(scaffolding_extenders));
    push_back_all(result, ExtractExtenders(loop_resolving_extenders));

    INFO("Using " << pe_libs << " paired-end " << support_.LibStr(pe_libs));
    INFO("Using " << scf_pe_libs << " paired-end scaffolding " << support_.LibStr(scf_pe_libs));
    INFO("Using " << single_read_libs << " single read " << support_.LibStr(single_read_libs));

    PrintExtenders(result);
    return result;
}

}
