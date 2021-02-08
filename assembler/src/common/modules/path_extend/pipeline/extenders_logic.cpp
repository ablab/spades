//
// Created by andrey on 14.11.16.
//

#include "extenders_logic.hpp"
#include "modules/path_extend/scaffolder2015/extension_chooser2015.hpp"


namespace path_extend {

using namespace debruijn_graph;
using namespace std;
using namespace omnigraph::de;

shared_ptr<ExtensionChooser> ExtendersGenerator::MakeLongReadsExtensionChooser(size_t lib_index,
                                                                               const GraphCoverageMap &read_paths_cov_map) const {
    auto long_reads_config = support_.GetLongReadsConfig(dataset_info_.reads[lib_index].type());

    if (dataset_info_.reads[lib_index].type() == io::LibraryType::TrustedContigs) {
        return make_shared<TrustedContigsExtensionChooser>(graph_, read_paths_cov_map,
                                                            long_reads_config.filtering,
                                                            long_reads_config.weight_priority,
                                                            long_reads_config.unique_edge_priority,
                                                            long_reads_config.min_significant_overlap,
                                                            params_.pset.extension_options.max_repeat_length,
                                                            params_.uneven_depth);
    }

    return make_shared<LongReadsExtensionChooser>(graph_, read_paths_cov_map,
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
        //TODO does max make sense here?
        resolvable_repeat_length_bound = std::max(resolvable_repeat_length_bound, lib.data().unmerged_read_length);
    }
    INFO("resolvable_repeat_length_bound set to " << resolvable_repeat_length_bound);
    bool investigate_short_loop = lib.is_contig_lib() || lib.is_long_read_lib() || support_.UseCoverageResolverForSingleReads(lib.type());

    auto long_read_ec = MakeLongReadsExtensionChooser(lib_index, read_paths_cov_map);
    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       used_unique_storage_,
                                       long_read_ec,
                                       resolvable_repeat_length_bound,
                                       investigate_short_loop, /* investigate short loops */
                                       support_.UseCoverageResolverForSingleReads(lib.type()));
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeLongEdgePEExtender(size_t lib_index,
                                                                      bool investigate_loops) const {
    const auto &clustered_indices = gp_.get<PairedInfoIndicesT<Graph>>("clustered_indices");

    const auto &lib = dataset_info_.reads[lib_index];
    auto paired_lib = MakeNewLib(graph_, lib, clustered_indices[lib_index]);
    //INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<WeightCounter> wc =
        make_shared<PathCoverWeightCounter>(graph_, paired_lib,
                                            params_.pset.normalize_weight,
                                            params_.pset.extension_options.single_threshold);
    auto opts = support_.GetExtensionOpts(paired_lib, params_.pset);
    shared_ptr<ExtensionChooser> extension =
        make_shared<LongEdgeExtensionChooser>(graph_, wc,
                                              opts.weight_threshold,
                                              opts.priority_coeff);

    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       used_unique_storage_,
                                       extension,
                                       paired_lib->GetISMax(),
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/,
                                       opts.weight_threshold);
}

shared_ptr<GapAnalyzer> ExtendersGenerator::MakeGapAnalyzer(double is_variation) const {
    const auto &pset = params_.pset;

    vector<shared_ptr<GapAnalyzer>> joiners;
    if (params_.pset.scaffolder_options.use_la_gap_joiner)
        joiners.push_back(std::make_shared<LAGapAnalyzer>(graph_, pset.scaffolder_options.min_overlap_length,
                                                        pset.scaffolder_options.flank_multiplication_coefficient,
                                                        pset.scaffolder_options.flank_addition_coefficient));


    joiners.push_back(std::make_shared<HammingGapAnalyzer>(graph_,
                                                         pset.scaffolder_options.min_gap_score,
                                                         pset.scaffolder_options.short_overlap,
                                                         (int) pset.scaffolder_options.basic_overlap_coeff
                                                             * dataset_info_.RL));

    //todo introduce explicit must_overlap_coeff and rename max_can_overlap -> can_overlap_coeff
    return std::make_shared<CompositeGapAnalyzer>(graph_,
                                                joiners,
                                                size_t(math::round(pset.scaffolder_options.max_can_overlap
                                                           * is_variation)), /* may overlap threshold */
                                                int(math::round(-pset.scaffolder_options.var_coeff * is_variation)), /* must overlap threshold */
                                                pset.scaffolder_options.artificial_gap);

}

shared_ptr<PathExtender> ExtendersGenerator::MakeScaffoldingExtender(size_t lib_index) const {

    const auto &lib = dataset_info_.reads[lib_index];
    const auto &pset = params_.pset;
    const auto &scaffolding_indices = gp_.get<PairedInfoIndicesT<Graph>>("scaffolding_indices");
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(graph_, lib, scaffolding_indices[lib_index]);

    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(graph_, paired_lib);

    auto scaff_chooser = std::make_shared<ScaffoldingExtensionChooser>(graph_, counter,
                                                                       pset.scaffolder_options.cl_threshold,
                                                                       pset.scaffolder_options.var_coeff);

    return make_shared<ScaffoldingPathExtender>(gp_, cover_map_,
                                                used_unique_storage_, scaff_chooser,
                                                MakeGapAnalyzer(paired_lib->GetIsVar()),
                                                paired_lib->GetISMax(),
                                                false, /* investigate short loops */
                                                params_.avoid_rc_connections);
}

shared_ptr<PathExtender> ExtendersGenerator::MakeRNAScaffoldingExtender(size_t lib_index) const {

    const auto &lib = dataset_info_.reads[lib_index];
    const auto &pset = params_.pset;
    const auto &paired_indices = gp_.get<UnclusteredPairedInfoIndicesT<Graph>>();
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(graph_, lib, paired_indices[lib_index]);

    shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(graph_, paired_lib);

    auto scaff_chooser = std::make_shared<ScaffoldingExtensionChooser>(graph_,
                                                                       counter,
                                                                       pset.scaffolder_options.cutoff,
                                                                       pset.scaffolder_options.var_coeff,
                                                                       pset.scaffolder_options.rel_cov_cutoff);
    auto scaff_chooser2 = std::make_shared<ScaffoldingExtensionChooser>(graph_,
                                                                        counter,
                                                                        pset.scaffolder_options.hard_cutoff,
                                                                        pset.scaffolder_options.var_coeff,
                                                                        pset.scaffolder_options.rel_cov_cutoff);


    VERIFY(pset.scaffolder_options.min_overlap_for_rna_scaffolding.is_initialized());
    return make_shared<RNAScaffoldingPathExtender>(gp_, cover_map_,
                                                   used_unique_storage_,
                                                   scaff_chooser,
                                                   scaff_chooser2,
                                                   MakeGapAnalyzer(paired_lib->GetIsVar()),
                                                   paired_lib->GetISMax(),
                                                   false  /* investigate short loops */,
                                                   *pset.scaffolder_options.min_overlap_for_rna_scaffolding);
}

std::shared_ptr<ExtensionChooser> ExtendersGenerator::MakeLongReadsRNAExtensionChooser(size_t lib_index,
                                                                                  const GraphCoverageMap &read_paths_cov_map) const {
    auto long_reads_config = support_.GetLongReadsConfig(dataset_info_.reads[lib_index].type());
    INFO("Creating long read rna chooser")
    return std::make_shared<LongReadsRNAExtensionChooser>(graph_, read_paths_cov_map,
                                                     long_reads_config.filtering,
                                                     long_reads_config.min_significant_overlap);
}


std::shared_ptr<SimpleExtender> ExtendersGenerator::MakeLongReadsRNAExtender(size_t lib_index,
                                                                        const GraphCoverageMap& read_paths_cov_map) const {
    const auto& lib = dataset_info_.reads[lib_index];
    //TODO params
    size_t resolvable_repeat_length_bound = 10000ul;
    if (!dataset_info_.reads[lib_index].is_contig_lib()) {
        resolvable_repeat_length_bound = std::max(resolvable_repeat_length_bound, lib.data().unmerged_read_length);
    }
    INFO("resolvable_repeat_length_bound set to " << resolvable_repeat_length_bound);
    bool investigate_short_loop = false;

    auto long_read_ec = MakeLongReadsRNAExtensionChooser(lib_index, read_paths_cov_map);
    INFO("Creating long read rna extender")
    return std::make_shared<MultiExtender>(gp_, cover_map_,
                                           used_unique_storage_,
                                           long_read_ec,
                                           resolvable_repeat_length_bound,
                                           investigate_short_loop, /* investigate short loops */
                                           false,
                                           0.0);
}

shared_ptr<PathExtender> ExtendersGenerator::MakeMatePairScaffoldingExtender(size_t lib_index,
                                                                             const ScaffoldingUniqueEdgeStorage &storage) const {

    const auto &lib = dataset_info_.reads[lib_index];
    const auto &pset = params_.pset;
    const auto &paired_indices = gp_.get<UnclusteredPairedInfoIndicesT<Graph>>();
    const auto &clustered_indices = gp_.get<PairedInfoIndicesT<Graph>>("clustered_indices");

    shared_ptr<PairedInfoLibrary> paired_lib;
    INFO("Creating Scaffolding 2015 extender for lib #" << lib_index);

    //FIXME: DimaA
    if (paired_indices[lib_index].size() > clustered_indices[lib_index].size()) {
        INFO("Paired unclustered indices not empty, using them");
        paired_lib = MakeNewLib(graph_, lib, paired_indices[lib_index]);
    } else if (clustered_indices[lib_index].size()) {
        INFO("clustered indices not empty, using them");
        paired_lib = MakeNewLib(graph_, lib, clustered_indices[lib_index]);
    } else {
        ERROR("All paired indices are empty!");
    }

    //TODO::was copypasted from MakeScaffoldingExtender, refactor 2015 extension chooser
    DEBUG("creating extchooser");
    shared_ptr<ConnectionCondition>
        condition = make_shared<PairedLibConnectionCondition>(graph_, paired_lib, lib_index, 0);
    auto scaff_chooser = std::make_shared<ExtensionChooser2015>(graph_,
                                                                nullptr,
                                                                condition,
                                                                storage,
                                                                pset.scaffolder_options.cl_threshold,
                                                                pset.scaffolder_options.var_coeff,
                                                                pset.scaffolding2015.relative_weight_cutoff,
                                                                graph_.size()
                                                                    <= params_.pset.scaffolding2015.graph_connectivity_max_edges);

    return make_shared<ScaffoldingPathExtender>(gp_, cover_map_,
                                                used_unique_storage_,
                                                scaff_chooser,
                                                MakeGapAnalyzer(paired_lib->GetIsVar()),
                                                paired_lib->GetISMax(),
                                                false, /* investigate short loops */
                                                params_.avoid_rc_connections,
                                                false /* jump only from tips */);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeCoordCoverageExtender(size_t lib_index) const {
    const auto& lib = dataset_info_.reads[lib_index];
    const auto &clustered_indices = gp_.get<PairedInfoIndicesT<Graph>>("clustered_indices");
    auto paired_lib = MakeNewLib(graph_, lib, clustered_indices[lib_index]);

    auto provider = make_shared<CoverageAwareIdealInfoProvider>(graph_, paired_lib, lib.data().unmerged_read_length);

    auto meta_wc = make_shared<PathCoverWeightCounter>(graph_, paired_lib,
                                                       params_.pset.normalize_weight,
                                                       params_.pset.extension_options.single_threshold,
                                                       provider);

    auto permissive_pi_chooser = make_shared<IdealBasedExtensionChooser>(graph_,
                                                                         meta_wc,
                                                                         params_.pset.extension_options.weight_threshold,
                                                                         params_.pset.extension_options.priority_coeff);

    auto coord_cov_chooser = make_shared<CoordinatedCoverageExtensionChooser>(graph_, *provider,
                                                                              params_.pset.coordinated_coverage.max_edge_length_in_repeat,
                                                                              params_.pset.coordinated_coverage.delta,
                                                                              params_.pset.coordinated_coverage.min_path_len);

    auto chooser = make_shared<JointExtensionChooser>(graph_, permissive_pi_chooser, coord_cov_chooser);

    return make_shared<SimpleExtender>(gp_, cover_map_, used_unique_storage_, chooser,
                                       -1ul /* insert size is needed only for loop detection, which is not needed in this case */,
                                       false, /* investigate short loops */
                                       false /*use short loop coverage resolver*/,
                                       params_.pset.extension_options.weight_threshold);
}

shared_ptr<SimpleExtender> ExtendersGenerator::MakeRNAExtender(size_t lib_index, bool investigate_loops) const {

    const auto &lib = dataset_info_.reads[lib_index];
    const auto &clustered_indices = gp_.get<PairedInfoIndicesT<Graph>>("clustered_indices");
    auto paired_lib = MakeNewLib(graph_, lib, clustered_indices[lib_index]);
//    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    auto cip = make_shared<CoverageAwareIdealInfoProvider>(graph_, paired_lib, lib.data().unmerged_read_length);
    shared_ptr<WeightCounter> wc =
        make_shared<PathCoverWeightCounter>(graph_, paired_lib, params_.pset.normalize_weight,
                                            params_.pset.extension_options.single_threshold,
                                            cip);

    auto opts = support_.GetExtensionOpts(paired_lib, params_.pset);
    shared_ptr<RNAExtensionChooser> extension =
        make_shared<RNAExtensionChooser>(graph_, wc,
                                         opts.weight_threshold,
                                         opts.priority_coeff);

    return make_shared<MultiExtender>(gp_, cover_map_,
                                      used_unique_storage_,
                                      extension,
                                      paired_lib->GetISMax(),
                                      investigate_loops,
                                      false /*use short loop coverage resolver*/,
                                      opts.weight_threshold);
}


shared_ptr<SimpleExtender> ExtendersGenerator::MakeSimpleCoverageExtender(size_t lib_index) const {

    const auto &ss_coverage = gp_.get<SSCoverageContainer>();
    auto extension =
        make_shared<SimpleCoverageExtensionChooser>(ss_coverage[lib_index], graph_,
                                                    params_.pset.simple_coverage_resolver.coverage_margin,
                                                    params_.pset.simple_coverage_resolver.max_coverage_variation,
                                                    params_.pset.simple_coverage_resolver.min_upper_coverage);

    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       used_unique_storage_,
                                       extension,
                                       1000, /*insert size for cycle detection*/
                                       false /*investigate short loops*/,
                                       false /*use short loop coverage resolver*/);
}


shared_ptr<SimpleExtender> ExtendersGenerator::MakePEExtender(size_t lib_index, bool investigate_loops) const {
    const auto &lib = dataset_info_.reads[lib_index];
    const auto &clustered_indices = gp_.get<PairedInfoIndicesT<Graph>>("clustered_indices");
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(graph_, lib, clustered_indices[lib_index]);
    VERIFY_MSG(!paired_lib->IsMp(), "Tried to create PE extender for MP library");
    auto opts = params_.pset.extension_options;
//    INFO("Threshold for lib #" << lib_index << ": " << paired_lib->GetSingleThreshold());

    shared_ptr<CoverageAwareIdealInfoProvider> iip = nullptr;
    if (params_.uneven_depth || params_.mode == config::pipeline_type::moleculo) {
        iip = make_shared<CoverageAwareIdealInfoProvider>(graph_, paired_lib, lib.data().unmerged_read_length);
    } else {
        double lib_cov = support_.EstimateLibCoverage(lib_index);
        INFO("Estimated coverage of library #" << lib_index << " is " << lib_cov);
        iip = make_shared<GlobalCoverageAwareIdealInfoProvider>(graph_, paired_lib, lib.data().unmerged_read_length, lib_cov);
    }

    auto wc = make_shared<PathCoverWeightCounter>(graph_, paired_lib, params_.pset.normalize_weight,
                                                  params_.pset.extension_options.single_threshold,
                                                  iip);

    auto extension_chooser = make_shared<SimpleExtensionChooser>(graph_, wc,
                                                         opts.weight_threshold,
                                                         opts.priority_coeff);
    INFO ("Creating extender; library index size: " << extension_chooser->wc()->PairedLibrary().size());

    return make_shared<SimpleExtender>(gp_, cover_map_,
                                       used_unique_storage_,
                                       extension_chooser,
                                       paired_lib->GetISMax(),
                                       investigate_loops,
                                       false /*use short loop coverage resolver*/,
                                       opts.weight_threshold);
}

//FIXME do we need ExtenderTriplets story here?
//FIXME integrate with MakeBasicExtenders
Extenders ExtendersGenerator::MakePEExtenders() const {
    Extenders result;
    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info_.reads[lib_index];
        if (support_.IsForPEExtender(lib)) {
            result.push_back(MakePEExtender(lib_index, false));
        }
    }
    return result;
}

void ExtendersGenerator::PrintExtenders(const Extenders &extenders) const {
    DEBUG("Extenders in vector:");
    for (const auto& extender : extenders) {
        //TODO: use polymorphism instead of RTTI
        auto ext_ptr = extender.get();
        DEBUG("Extender #i" << typeid(*ext_ptr).name());
        if (utils::instanceof<SimpleExtender>(ext_ptr)) {
            auto ec = ((SimpleExtender *) ext_ptr)->GetExtensionChooser();
            auto ec_ptr = ec.get();
            DEBUG("    Extender #i" << typeid(*ec_ptr).name());
        }
        else if (utils::instanceof<ScaffoldingPathExtender>(ext_ptr)) {
            auto ec = ((ScaffoldingPathExtender *) ext_ptr)->GetExtensionChooser();
            auto ec_ptr = ec.get();
            DEBUG("    Extender #i" << typeid(*ec_ptr).name());
        }
    }
}

Extenders ExtendersGenerator::MakeMPExtenders() const {
    Extenders extenders = MakeMPExtenders(unique_data_.main_unique_storage_);
    INFO("Using " << extenders.size() << " mate-pair " << support_.LibStr(extenders.size()));

    for (const auto& unique_storage : unique_data_.unique_storages_) {
        utils::push_back_all(extenders, MakeMPExtenders(unique_storage));
    }
    return extenders;
}

Extenders ExtendersGenerator::MakeMPExtenders(const ScaffoldingUniqueEdgeStorage &storage) const {
    ExtenderTriplets result;

    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info_.reads[lib_index];

        if (lib.is_mate_pair()) {
            result.emplace_back(lib.type(), lib_index,
                                MakeMatePairScaffoldingExtender(lib_index, storage));
        }
    }
    std::stable_sort(result.begin(), result.end());

    return ExtractExtenders(result);
}

Extenders ExtendersGenerator::MakePBScaffoldingExtenders() const {
    const auto &pset = params_.pset;
    ExtenderTriplets result;

    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); lib_index++) {
        if (support_.IsForSingleReadScaffolder(dataset_info_.reads[lib_index])) {
            INFO("Creating scaffolding extender for lib " << lib_index);
            shared_ptr<ConnectionCondition> condition = make_shared<LongReadsLibConnectionCondition>(graph_,
                                                                                                     lib_index, 2,
                                                                                                     unique_data_.long_reads_cov_map_[lib_index]);
            auto scaff_chooser = std::make_shared<ExtensionChooser2015>(graph_,
                                                                        nullptr,
                                                                        condition,
                                                                        unique_data_.unique_pb_storage_,
                                                                        pset.scaffolder_options.cl_threshold,
                                                                        pset.scaffolder_options.var_coeff,
                                                                        pset.scaffolding2015.relative_weight_cutoff);

            result.emplace_back(dataset_info_.reads[lib_index].type(),
                                lib_index,
                                //FIXME are utilized constants reasonable?
                                make_shared<ScaffoldingPathExtender>(gp_, cover_map_,
                                                                     used_unique_storage_,
                                                                     scaff_chooser,
                                                                     MakeGapAnalyzer(1000), /* "IS variation" */
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

Extenders ExtendersGenerator::MakeBasicExtenders() const {
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
        if (support_.IsForSingleReadExtender(lib) && ! cfg::get().pd) {
            if (!config::PipelineHelper::IsPlasmidPipeline(params_.mode)) {

                if (pset.multi_path_extend) {
                    basic_extenders.emplace_back(lib.type(), lib_index, MakeLongReadsRNAExtender(lib_index,
                                                                                             unique_data_.long_reads_cov_map_[lib_index]));
                    INFO("Created for lib #" << lib_index);
                } else {
                    basic_extenders.emplace_back(lib.type(), lib_index,
                                                 MakeLongReadsExtender(lib_index,
                                                                       unique_data_.long_reads_cov_map_[lib_index]));
                }
                ++single_read_libs;
            } else {
                if (lib.type() != io::LibraryType::PathExtendContigs)
                    WARN("Long reads and contig libraries are currently not supported within plasmid pipeline");
            }
        }
        if (support_.IsForPEExtender(lib)) {
            ++pe_libs;
            if (IsOldPEEnabled(pset.sm)) {
                if (params_.mode == config::pipeline_type::moleculo) {
                    basic_extenders.emplace_back(lib.type(), lib_index, MakeLongEdgePEExtender(lib_index, false));
                } else if (pset.multi_path_extend) {
                    basic_extenders.emplace_back(lib.type(), lib_index, MakePEExtender(lib_index, false));
                    if (pset.simple_coverage_resolver.enabled)
                        basic_extenders.emplace_back(lib.type(), lib_index, MakeSimpleCoverageExtender(lib_index));
                    basic_extenders.emplace_back(lib.type(), lib_index, MakeRNAExtender(lib_index, false));
                } else {
                    basic_extenders.emplace_back(lib.type(), lib_index, MakePEExtender(lib_index, false));
                }
            } else if (pset.sm == scaffolding_mode::sm_2015) {
                basic_extenders.emplace_back(lib.type(), lib_index,
                                             MakeMatePairScaffoldingExtender(lib_index,
                                                                             unique_data_.main_unique_storage_));
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
//Do we really need this in plasmid modes?
                scaffolding_extenders.emplace_back(lib.type(), lib_index, MakeScaffoldingExtender(lib_index));
                if (pset.sm == scaffolding_mode::sm_combined) {
                    scaffolding_extenders.emplace_back(lib.type(), lib_index,
                                                       MakeMatePairScaffoldingExtender(lib_index,
                                                                                       unique_data_.main_unique_storage_));
                }
            }
        }
    }

    std::stable_sort(basic_extenders.begin(), basic_extenders.end());
    std::stable_sort(scaffolding_extenders.begin(), scaffolding_extenders.end());
    std::stable_sort(loop_resolving_extenders.begin(), loop_resolving_extenders.end());

    Extenders result;
    utils::push_back_all(result, ExtractExtenders(basic_extenders));
    utils::push_back_all(result, ExtractExtenders(scaffolding_extenders));
    utils::push_back_all(result, ExtractExtenders(loop_resolving_extenders));

    INFO("Using " << pe_libs << " paired-end " << support_.LibStr(pe_libs));
    INFO("Using " << scf_pe_libs << " paired-end scaffolding " << support_.LibStr(scf_pe_libs));
    INFO("Using " << single_read_libs << " single read " << support_.LibStr(single_read_libs));

    PrintExtenders(result);
    return result;
}

}
