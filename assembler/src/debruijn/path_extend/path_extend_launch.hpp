//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
#include "pe_io.hpp"
#include "path_visualizer.hpp"
#include "path_validator.hpp"
#include "loop_traverser.h"
#include "single_threshold_finder.hpp"

//#include "gap_jumper.hpp"

#include "../include/omni/edges_position_handler.hpp"

namespace path_extend {

using namespace debruijn_graph;


/*void resolve_repeats_pe(size_t k, conj_graph_pack& gp, PairedInfoLibraries& libs, PairedInfoLibraries& scafolding_libs,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    INFO("Path extend repeat resolving tool started");
    params = p;

    make_dir(output_dir);
    std::string etcDir = output_dir + "path_extend/";
    if (cfg::get().make_saves || cfg::get().run_mode) {
        make_dir(etcDir);
    }

    if (!cfg::get().run_mode) {
        cfg::get().pe_params.debug_output = false;
        cfg::get().pe_params.output.DisableAll();
        cfg::get().pe_params.viz.DisableAll();
    }

    DEBUG("Using " << cfg::get().pe_params.param_set.metric << " metric");
    DEBUG("Scaffolder is " << (cfg::get().pe_params.param_set.scaffolder_options.on ? "on" : "off"));

    PathInfoWriter path_writer;
	PathVisualizer visualizer(k);
	ContigWriter writer(gp, k);

    if (cfg::get().pe_params.debug_output) {
        writer.writeEdges(etcDir + "before_resolve.fasta");
        visualizer.writeGraphSimple(gp, etcDir + "before_resolve.dot", "before_resolve");
    }


    DEBUG("Initializing weight counters");
	WeightCounter * wc = 0;
    WeightCounter * scaf_wc = 0;

    scaf_wc = new ReadCountWeightCounter(gp.g, scafolding_libs, cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold);
	if (cfg::get().pe_params.param_set.metric == "read_count") {
	    wc = new ReadCountWeightCounter(gp.g, libs, cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold);
	}
	else if (cfg::get().pe_params.param_set.metric == "path_cover") {
	    wc = new PathCoverWeightCounter(gp.g, libs, cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold, *cfg::get().pe_params.param_set.extension_options.select_options.single_threshold);
	} else {
	    WARN("Unknown metric");
	    return;
	}

	wc->setNormalizeWeight(cfg::get().pe_params.param_set.normalize_weight);
	wc->setNormalizeWightByCoverage(cfg::get().pe_params.param_set.normalize_by_coverage);

	DEBUG("Initializing resolver");
	ExtensionChooser * seedEC;
	if (cfg::get().pe_params.param_set.seed_selection.check_trusted) {
	    ExtensionChooser * trivialEC = new TrivialExtensionChooser(gp.g);
	    ExtensionChooser * simpleEC = new SimpleExtensionChooser(gp.g, wc, 1);
	    simpleEC->setExcludeTrivial(false);
        simpleEC->setExcludeTrivialWithBulges(false);

	    seedEC = new JointExtensionChooser(gp.g, trivialEC, simpleEC);
	} else {
	    seedEC = new TrivialExtensionChooser(gp.g);
	}

	CoveringPathExtender * seedPE = new SimplePathExtender(gp.g, 1, seedEC);
	seedPE->setInvestigateShortLoops(false);
	PathExtendResolver resolver(gp.g, k);

	INFO("Finding seeds");
	auto seeds = resolver.makeSeeds(*seedPE, cfg::get().pe_params.param_set.seed_selection.start_egde_coverage);
	INFO("Found " << seeds.size() << " seeds");


	if (cfg::get().pe_params.output.write_seeds) {
	    writer.writePaths(seeds, etcDir + "seeds.fasta");
	}
    if (cfg::get().pe_params.viz.print_seeds) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "seeds.dot", "Seeds", seeds);
        path_writer.writePaths(seeds, etcDir + "seeds.data");
    }

    if ( cfg::get().pe_params.param_set.seed_selection.min_coverage > 0) {
        DEBUG("Filtering seeds by coverage " << cfg::get().pe_params.param_set.seed_selection.min_coverage);
        CoveragePathFilter coverageFilter(gp.g, cfg::get().pe_params.param_set.seed_selection.min_coverage);
        coverageFilter.filter(seeds);
    }

    SimpleExtensionChooser * extensionChooser = new SimpleExtensionChooser(gp.g, wc, cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);
    SimplePathExtender * usualPE = new SimplePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, extensionChooser);

    ScaffoldingExtensionChooser * scafExtensionChooser = new ScaffoldingExtensionChooser(gp.g, scaf_wc, cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);
    GapJoiner * gapJoiner = new HammingGapJoiner(gp.g, cfg::get().pe_params.param_set.scaffolder_options.min_gap_score,
            (int) (cfg::get().pe_params.param_set.scaffolder_options.max_must_overlap * gp.g.k()),
            (int) (cfg::get().pe_params.param_set.scaffolder_options.max_can_overlap * gp.g.k()),
            cfg::get().pe_params.param_set.scaffolder_options.short_overlap);
	ScaffoldingPathExtender * scafPE = new ScaffoldingPathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, scafExtensionChooser, gapJoiner);

	CoveringPathExtender * mainPE = new CompositePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, usualPE, scafPE);

	seeds.SortByLength();
	INFO("Extending seeds");

	auto paths = resolver.extendSeeds(seeds, *mainPE);
    if (cfg::get().pe_params.output.write_overlaped_paths) {
        writer.writePaths(paths, etcDir + "overlaped_paths.fasta");
    }
    if (cfg::get().pe_params.viz.print_overlaped_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "overlaped_paths.dot", "Overlaped_paths", paths);
        path_writer.writePaths(paths, etcDir + "overlaped_paths.data");
    }

//  mainPE->GetCoverageMap().PrintMulticovered();

	if (cfg::get().pe_params.param_set.filter_options.remove_overlaps) {
	    INFO("Removing overlaps");
	    resolver.removeOverlaps(paths, mainPE->GetCoverageMap());
	}

	//mainPE->GetCoverageMap().PrintUncovered();
	//paths.CheckSymmetry();

	DEBUG("Adding uncovered edges");
	resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());

	paths.SortByLength();

	INFO("Found " << paths.size() << " contigs");
	writer.writePaths(paths, output_dir + contigs_name);
    if (cfg::get().pe_params.output.write_paths) {
        writer.writePathEdges(paths, output_dir + "paths.dat");
    }

    if (cfg::get().pe_params.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "paths.dot", "Final_paths", paths);
        path_writer.writePaths(paths, etcDir + "paths.data");
    }

	INFO("Path extend repeat resolving tool finished");
}


void resolve_repeats_pe_wcontigs(size_t k, conj_graph_pack& gp, PairedInfoLibraries& libs, PairedInfoLibraries& scafolding_libs,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    INFO("Path extend repeat resolving tool started");
    params = p;

    make_dir(output_dir);
    std::string etcDir = output_dir + "path_extend/";
    if (cfg::get().make_saves || cfg::get().run_mode) {
        make_dir(etcDir);
    }

    if (!cfg::get().run_mode) {
        cfg::get().pe_params.debug_output = false;
        cfg::get().pe_params.output.DisableAll();
        cfg::get().pe_params.viz.DisableAll();
    }

    DEBUG("Using " << cfg::get().pe_params.param_set.metric << " metric");
    DEBUG("Scaffolder is " << (cfg::get().pe_params.param_set.scaffolder_options.on ? "on" : "off"));

    PathInfoWriter path_writer;
    PathVisualizer visualizer(k);
    ContigWriter writer(gp, k);

    if (cfg::get().pe_params.debug_output) {
        writer.writeEdges(etcDir + "before_resolve.fasta");
        visualizer.writeGraphSimple(gp, etcDir + "before_resolve.dot", "before_resolve");
    }

    PathContainer supportingContigs;
    INFO("PARAMS FOR ADDITIONAL_CONTIGS" << cfg::get().pe_params.additional_contigs);
    if (fileExists(cfg::get().pe_params.additional_contigs)) {
        INFO("Reading additional contigs !!!!!!!" << cfg::get().pe_params.additional_contigs);
        supportingContigs = CreatePathsFromContigs(gp, cfg::get().pe_params.additional_contigs);
        writer.writePaths(supportingContigs, etcDir + "supporting_paths.fasta");
    }

    DEBUG("Initializing weight counters");
    WeightCounter * wc = 0;
    WeightCounter * scaf_wc = 0;

    scaf_wc = new ReadCountWeightCounter(gp.g, scafolding_libs, cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold);
    if (cfg::get().pe_params.param_set.metric == "read_count") {
        wc = new ReadCountWeightCounter(gp.g, libs, cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold);
    }
    else if (cfg::get().pe_params.param_set.metric == "path_cover") {
        wc = new PathCoverWeightCounter(gp.g, libs, cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold, *cfg::get().pe_params.param_set.extension_options.select_options.single_threshold);
    } else {
        WARN("Unknown metric");
        return;
    }

    wc->setNormalizeWeight(cfg::get().pe_params.param_set.normalize_weight);
    wc->setNormalizeWightByCoverage(cfg::get().pe_params.param_set.normalize_by_coverage);

    DEBUG("Initializing resolver");
    ExtensionChooser * seedEC;
    if (cfg::get().pe_params.param_set.seed_selection.check_trusted) {
        ExtensionChooser * trivialEC = new TrivialExtensionChooser(gp.g);
        ExtensionChooser * simpleEC = new SimpleExtensionChooser(gp.g, wc, 1);
        simpleEC->setExcludeTrivial(false);
        simpleEC->setExcludeTrivialWithBulges(false);

        seedEC = *//*new TrivialExtensionChooserWithPI(gp.g, wc);*//*new JointExtensionChooser(gp.g, trivialEC, simpleEC);
    } else {
        seedEC = new TrivialExtensionChooser(gp.g);
    }

    LoopDetectingPathExtender * seedPE = new SimplePathExtender(gp.g, 1, seedEC);
    seedPE->setInvestigateShortLoops(false);
    PathExtendResolver resolver(gp.g, k);

    INFO("Finding seeds");
    auto seeds = resolver.makeSeeds(*seedPE, cfg::get().pe_params.param_set.seed_selection.start_egde_coverage);//ReturnEdges(0);
    INFO("Found " << seeds.size() << " seeds");


    if (cfg::get().pe_params.output.write_seeds) {
        writer.writePaths(seeds, etcDir + "seeds.fasta");
    }
    if (cfg::get().pe_params.viz.print_seeds) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "seeds.dot", "Seeds", seeds);
        path_writer.writePaths(seeds, etcDir + "seeds.data");
    }

    if ( cfg::get().pe_params.param_set.seed_selection.min_coverage > 0) {
        DEBUG("Filtering seeds by coverage " << cfg::get().pe_params.param_set.seed_selection.min_coverage);
        CoveragePathFilter coverageFilter(gp.g, cfg::get().pe_params.param_set.seed_selection.min_coverage);
        coverageFilter.filter(seeds);
    }
    SimpleExtensionChooser * extensionChooser = new SimpleExtensionChooser(gp.g, wc, cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);
    SimplePathExtender * usualPE = new SimplePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, extensionChooser);

    ScaffoldingExtensionChooser * scafExtensionChooser = new ScaffoldingExtensionChooser(gp.g, scaf_wc, cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);
    GapJoiner * gapJoiner = new HammingGapJoiner(gp.g, cfg::get().pe_params.param_set.scaffolder_options.min_gap_score,
            (int) (cfg::get().pe_params.param_set.scaffolder_options.max_must_overlap * gp.g.k()),
            (int) (cfg::get().pe_params.param_set.scaffolder_options.max_can_overlap * gp.g.k()),
            cfg::get().pe_params.param_set.scaffolder_options.short_overlap);
    ScaffoldingPathExtender * scafPE = new ScaffoldingPathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, scafExtensionChooser, gapJoiner);

    ExtensionChooser * pdEC = new PathsDrivenExtensionChooser(gp.g, supportingContigs);
    SimplePathExtender * pdPE = new SimplePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, pdEC);
    CoveringPathExtender * mainPE = new CompositePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, usualPE, scafPE, pdPE); //only pdPE


    seeds.SortByLength();
    INFO("Extending seeds");

    auto paths = resolver.extendSeeds(seeds, *mainPE);
    if (params.output.write_overlaped_paths) {
        writer.writePaths(paths, etcDir + "overlaped_paths.fasta");
    }
    if (cfg::get().pe_params.viz.print_overlaped_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "overlaped_paths.dot", "Overlaped_paths", paths);
        path_writer.writePaths(paths, etcDir + "overlaped_paths.data");
    }

//    mainPE->GetCoverageMap().PrintMulticovered();

    if (cfg::get().pe_params.param_set.filter_options.remove_overlaps) {
        INFO("Removing overlaps");
        resolver.removeOverlaps(paths, mainPE->GetCoverageMap());
    }

//    mainPE->GetCoverageMap().PrintMulticovered();
//  mainPE->GetCoverageMap().PrintUncovered();
    paths.CheckSymmetry();

    DEBUG("Adding uncovered edges");
    resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());
    INFO("loop Traverser");
    LoopTraverser loopTraverser(gp.g, paths, mainPE->GetCoverageMap());
    loopTraverser.TraverseAllLoops();

    paths.SortByLength();

    INFO("Found " << paths.size() << " contigs");
    writer.writePaths(paths, output_dir + contigs_name);
    if (cfg::get().pe_params.output.write_paths) {
        writer.writePathEdges(paths, output_dir + "paths.dat");
    }

    if (params.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "paths.dot", "Final_paths", paths);
        path_writer.writePaths(paths, etcDir + "paths.data");
    }

    INFO("Path extend repeat resolving tool finished");

    if (supportingContigs.size() == paths.size()) {
        INFO("Comparing mapped shit")
        for (size_t i = 0; i < supportingContigs.size(); ++i) {
            if (supportingContigs[i] != paths[i]) {
                INFO("Non-equal paths. Before:");
                paths[i].Print();
                INFO("After")
                supportingContigs[i].Print();

            }
        }
    }
}*/




void resolve_repeats_pe_mp(size_t k, conj_graph_pack& gp,
		PairedInfoLibraries& libs, PairedInfoLibraries& scafolding_libs,
		PairedInfoLibraries& mate_pair_libs,
		PairedInfoLibraries& mate_pair_scafolding_libs,
		const std::string& output_dir, const std::string& contigs_name) {
	make_dir(output_dir);
	std::string etcDir = output_dir + "path_extend/";
	make_dir(etcDir);
	if (!cfg::get().run_mode) {
		params.output.DisableAll();
		params.viz.DisableAll();
	}
	ofstream out("./scaffolder.log", ostream::out);
	out.close();

	INFO("Path extend repeat resolving tool started");
	INFO("Using " << cfg::get().pe_params.param_set.metric << " metric");
	INFO(
			"Scaffolder is "
					<< (cfg::get().pe_params.param_set.scaffolder_options.on ?
							"on" : "off"));

	PathInfoWriter path_writer;
	PathVisualizer visualizer(k);
	ContigWriter writer(gp, k, libs[0]->insert_size_);
	writer.writeEdges(etcDir + "before_resolve.fasta");

	INFO("Initializing weight counters");
	WeightCounter * wc = 0;
	WeightCounter * mp_wc = 0;
	WeightCounter * scaf_wc = 0;

	scaf_wc = new ReadCountWeightCounter(gp.g, scafolding_libs);
	mp_wc =
			new PathCoverWeightCounter(gp.g, mate_pair_libs,
					cfg::get().pe_params.param_set.mate_pair_options.select_options.weight_threshold,
					*cfg::get().pe_params.param_set.mate_pair_options.select_options.single_threshold);

	if (cfg::get().pe_params.param_set.metric == "read_count") {
		wc =
				new ReadCountWeightCounter(gp.g, libs,
						cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold);
	} else if (cfg::get().pe_params.param_set.metric == "path_cover") {
		wc =
				new PathCoverWeightCounter(gp.g, libs,
						cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold,
						*cfg::get().pe_params.param_set.extension_options.select_options.single_threshold);
	} else {
		WARN("Unknown metric");
		return;
	}

	wc->setNormalizeWeight(cfg::get().pe_params.param_set.normalize_weight);
	wc->setNormalizeWightByCoverage(cfg::get().pe_params.param_set.normalize_by_coverage);
	mp_wc->setNormalizeWeight(cfg::get().pe_params.param_set.normalize_weight);
	mp_wc->setNormalizeWightByCoverage(cfg::get().pe_params.param_set.normalize_by_coverage);

	/*INFO("Initializing resolver");
	ExtensionChooser * seedEC;
	if (cfg::get().pe_params.param_set.seed_selection.check_trusted) {
		seedEC = new TrivialExtensionChooserWithPI(gp.g, wc);
	} else {
		seedEC = new TrivialExtensionChooser(gp.g);
	}

	PathExtender * seedPE = new SimplePathExtender(gp.g, 1, seedEC);*/
	PathExtendResolver resolver(gp.g, k);
	auto seeds = resolver.makeSimpleSeeds();
	/*INFO("Finding seeds");
	auto seeds = resolver.makeSeeds(*seedPE,
			cfg::get().pe_params.param_set.seed_selection.start_egde_coverage);
	INFO("Found " << seeds.size() << " seeds");

	if (params.output.write_seeds) {
		writer.writePaths(seeds, etcDir + "seeds.fasta");
	}
	if (params.viz.print_seeds) {
		visualizer.writeGraphWithPathsSimple(gp, etcDir + "seeds.dot", "Seeds",
				seeds);
		path_writer.writePaths(seeds, etcDir + "seeds.data");
	}

	if (cfg::get().pe_params.param_set.seed_selection.min_coverage > 0) {
		INFO(
				"Filtering seeds by coverage "
						<< cfg::get().pe_params.param_set.seed_selection.min_coverage);
		CoveragePathFilter coverageFilter(gp.g,
				cfg::get().pe_params.param_set.seed_selection.min_coverage);
		coverageFilter.filter(seeds);
	}*/

	SimpleExtensionChooser * extensionChooser = new SimpleExtensionChooser(gp.g,
			wc,
			cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);
	SimplePathExtender * usualPE = new SimplePathExtender(gp.g,
			cfg::get().pe_params.param_set.loop_removal.max_loops, extensionChooser);
	ScaffoldingExtensionChooser * scafExtensionChooser =
			new ScaffoldingExtensionChooser(gp.g, scaf_wc,
					cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);
	GapJoiner * gapJoiner = new HammingGapJoiner(gp.g,
			cfg::get().pe_params.param_set.scaffolder_options.min_gap_score,
			(int) (cfg::get().pe_params.param_set.scaffolder_options.max_must_overlap
					* gp.g.k()),
			(int) (cfg::get().pe_params.param_set.scaffolder_options.max_can_overlap
					* gp.g.k()),
			cfg::get().pe_params.param_set.scaffolder_options.short_overlap);
	ScaffoldingPathExtender * scafPE = new ScaffoldingPathExtender(gp.g,
			cfg::get().pe_params.param_set.loop_removal.max_loops, scafExtensionChooser,
			gapJoiner);
	SimpleExtensionChooser * mpExtensionChooser = new SimpleExtensionChooser(
			gp.g, mp_wc,
			cfg::get().pe_params.param_set.mate_pair_options.select_options.priority_coeff);
	SimplePathExtender * mpPE = new SimplePathExtender(gp.g,
			cfg::get().pe_params.param_set.loop_removal.max_loops, mpExtensionChooser);
	CoveringPathExtender * mainPE = new CompositePathExtender(gp.g,
				cfg::get().pe_params.param_set.loop_removal.max_loops, usualPE, mpPE);
	if (cfg::get().pe_params.param_set.scaffolder_options.on){
		INFO("Scaffolding is on, use scaf lib	");
		mainPE = new CompositePathExtender(gp.g,
		cfg::get().pe_params.param_set.loop_removal.max_loops, usualPE, scafPE, mpPE);
	}

	seeds.SortByLength();
	INFO("Extending seeds");

	auto paths = resolver.extendSeeds(seeds, *mainPE);
	if (cfg::get().pe_params.output.write_overlaped_paths) {
		writer.writePaths(paths, etcDir + "overlaped_paths.fasta");
	}
	if (cfg::get().pe_params.viz.print_overlaped_paths) {
		visualizer.writeGraphWithPathsSimple(gp, etcDir + "overlaped_paths.dot",
				"Overlaped_paths", paths);
		path_writer.writePaths(paths, etcDir + "overlaped_paths.data");
	}

	if (cfg::get().pe_params.param_set.filter_options.remove_overlaps) {
		INFO("Removing overlaps");
		resolver.removeOverlaps(paths, mainPE->GetCoverageMap(), mate_pair_libs[0]->insert_size_ + 2 * mate_pair_libs[0]->is_variation_);
	}

	mainPE->GetCoverageMap().PrintUncovered();
	paths.CheckSymmetry();

	INFO("loop Traverser");
	LoopTraverser loopTraverser(gp.g, paths, mainPE->GetCoverageMap(), mainPE);
	loopTraverser.TraverseAllLoops();
	paths.SortByLength();

	INFO("Adding uncovered edges");
	resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());

	INFO("Found " << paths.size() << " contigs");
	writer.writePaths(paths, output_dir + contigs_name);
	if (cfg::get().pe_params.output.write_paths) {
		writer.writePathEdges(paths, output_dir + "paths.dat");
	}

	if (cfg::get().pe_params.viz.print_paths) {
		visualizer.writeGraphWithPathsSimple(gp, etcDir + "paths.dot",
				"Final_paths", paths);
		path_writer.writePaths(paths, etcDir + "paths.data");
	}
	INFO("Path extend repeat resolving tool finished");
}
void resolve_repeats_pe_many_libs(size_t k, conj_graph_pack& gp,
		vector<PairedInfoLibraries>& libes, vector<PairedInfoLibraries>& scafolding_libes,
		const std::string& output_dir, const std::string& contigs_name) {
	INFO("Path extend repeat resolving tool started MANY LIBS");
	make_dir(output_dir);
	std::string etcDir = output_dir + "path_extend/";
	make_dir(etcDir);

	if (!cfg::get().run_mode) {
		cfg::get().pe_params.output.DisableAll();
		cfg::get().pe_params.viz.DisableAll();
	}
	ofstream out("./scaffolder.log", ostream::out);
	out.close();
	size_t maxOverlapedLength = 0;
	size_t minInsertSize = -1;

	for (size_t i = 0; i < libes.size(); ++i) {
		for (size_t j = 0; j < libes[i].size(); ++j) {
			size_t overlap = libes[i][j]->insert_size_ + libes[i][j]->is_variation_;
			if (overlap > maxOverlapedLength) {
				maxOverlapedLength = overlap;
			}
			if (libes[i][j]->insert_size_ < minInsertSize) {
				minInsertSize = libes[i][j]->insert_size_;
			}
		}
	}
	INFO("max overlaped length " << maxOverlapedLength);

	INFO("Using " << cfg::get().pe_params.param_set.metric << " metric");
	INFO("Scaffolder is " << (cfg::get().pe_params.param_set.scaffolder_options.on ? "on" : "off"));

	PathInfoWriter path_writer;
	PathVisualizer visualizer(k);
	ContigWriter writer(gp, k, minInsertSize);
	if (cfg::get().pe_params.debug_output) {
		writer.writeEdges(etcDir + "before_resolve.fasta");
		visualizer.writeGraphSimple(gp, etcDir + "before_resolve.dot", "before_resolve");
	}
	INFO("Initializing weight counters");
	vector<WeightCounter*> wcs;
	vector<WeightCounter*> scaf_wcs;
	PathContainer supportingContigs;
	if (FileExists(cfg::get().pe_params.additional_contigs)) {
		INFO("Reading additional contigs !!!!!!!" << cfg::get().pe_params.additional_contigs);
		supportingContigs = CreatePathsFromContigs(gp, cfg::get().pe_params.additional_contigs);
		writer.writePaths(supportingContigs, etcDir + "supporting_paths.fasta");
	}


	if (cfg::get().pe_params.param_set.metric == "read_count") {
		for (size_t i = 0; i < libes.size(); ++i){
			wcs.push_back(new ReadCountWeightCounter(gp.g, libes[i],
					cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold));
		}
	} else if (cfg::get().pe_params.param_set.metric == "path_cover") {
		for (size_t i = 0; i < libes.size(); ++i){
			wcs.push_back(new PathCoverWeightCounter(gp.g, libes[i],
						cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold,
						*cfg::get().pe_params.param_set.extension_options.select_options.single_threshold));
		}
	} else {
		WARN("Unknown metric");
		return;
	}
	vector<SimplePathExtender *> usualPEs;
	for (size_t i = 0; i < wcs.size(); ++i){
		wcs[i]->setNormalizeWeight(cfg::get().pe_params.param_set.normalize_weight);
		wcs[i]->setNormalizeWightByCoverage(cfg::get().pe_params.param_set.normalize_by_coverage);
		double priory_coef = libes[i][0]->is_mate_pair_?cfg::get().pe_params.param_set.mate_pair_options.select_options.priority_coeff : cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff;
		SimpleExtensionChooser * extensionChooser = new SimpleExtensionChooser(gp.g, wcs[i], priory_coef);
		usualPEs.push_back( new SimplePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, extensionChooser));
	}
	GapJoiner * gapJoiner = new HammingGapJoiner(gp.g,
				cfg::get().pe_params.param_set.scaffolder_options.min_gap_score,
				(int) (cfg::get().pe_params.param_set.scaffolder_options.max_must_overlap
						* gp.g.k()),
				(int) (cfg::get().pe_params.param_set.scaffolder_options.max_can_overlap
						* gp.g.k()),
				cfg::get().pe_params.param_set.scaffolder_options.short_overlap);
	vector<ScaffoldingPathExtender*> scafPEs;
	for (size_t i = 0; i < scafolding_libes.size(); ++i){
		scaf_wcs.push_back(new ReadCountWeightCounter(gp.g, scafolding_libes[i]));
		ScaffoldingExtensionChooser * scafExtensionChooser =
					new ScaffoldingExtensionChooser(gp.g, scaf_wcs[i],
							scafolding_libes[i][0]->is_mate_pair_?cfg::get().pe_params.param_set.mate_pair_options.select_options.priority_coeff : cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);
		scafPEs.push_back(new ScaffoldingPathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, scafExtensionChooser, gapJoiner));
	}
	PathExtendResolver resolver(gp.g, k);
	auto seeds = resolver.makeSimpleSeeds();
	ExtensionChooser * pdEC = new PathsDrivenExtensionChooser(gp.g, supportingContigs);
	SimplePathExtender * pdPE = new SimplePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, pdEC);

	vector<PathExtender *> all_libs(usualPEs.begin(), usualPEs.end());
	all_libs.insert(all_libs.end(), scafPEs.begin(), scafPEs.end());
	all_libs.push_back(pdPE);
	CoveringPathExtender * mainPE = new CompositePathExtender(gp.g, cfg::get().pe_params.param_set.loop_removal.max_loops, all_libs);

	seeds.SortByLength();
	INFO("Extending seeds");

	auto paths = resolver.extendSeeds(seeds, *mainPE);
	if (cfg::get().pe_params.output.write_overlaped_paths) {
		writer.writePaths(paths, etcDir + "overlaped_paths.fasta");
	}
	if (cfg::get().pe_params.viz.print_overlaped_paths) {
		visualizer.writeGraphWithPathsSimple(gp, etcDir + "overlaped_paths.dot",
				"Overlaped_paths", paths);
		path_writer.writePaths(paths, etcDir + "overlaped_paths.data");
	}
	INFO("before overlapes");
    for (size_t i = 0; i < paths.size(); ++i){
        paths.Get(i)->Print();
    }

    resolver.removeOverlaps(paths, mainPE->GetCoverageMap(), maxOverlapedLength);
	mainPE->GetCoverageMap().PrintUncovered();
	paths.CheckSymmetry();

    resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());
	paths.SortByLength();
    //writer.writePaths(paths, output_dir + contigs_name + "before_tr_loops.fasta");
	INFO("after overlapes");
    for (size_t i = 0; i < paths.size(); ++i){
        paths.Get(i)->Print();
    }
    paths.SortByLength();
    writer.writePaths(paths, output_dir + contigs_name);
	INFO("loop Traverser");
    LoopTraverser loopTraverser(gp.g, paths, mainPE->GetCoverageMap(), mainPE);
	loopTraverser.TraverseAllLoops();
	paths.SortByLength();
	INFO("Found " << paths.size() << " contigs");
	writer.writePaths(paths, output_dir + contigs_name.substr(0, contigs_name.rfind(".fasta")) + "_loop_tr.fasta");
	if (cfg::get().pe_params.output.write_paths) {
		writer.writePathEdges(paths, output_dir + "paths.dat");
	}

	if (cfg::get().pe_params.viz.print_paths) {
		visualizer.writeGraphWithPathsSimple(gp, etcDir + "paths.dot",
				"Final_paths", paths);
		path_writer.writePaths(paths, etcDir + "paths.data");
	}
	INFO("Path extend repeat resolving tool finished");
}

void resolve_repeats_pe(size_t k, conj_graph_pack& gp, PairedInfoIndexT<Graph>& paired_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(), paired_index));

    PairedInfoLibraries scaf_libs;
    vector<PairedInfoLibraries> libes;
    libes.push_back(libs);
    vector<PairedInfoLibraries> scafolding_libes;

    //resolve_repeats_pe_wcontigs(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
    resolve_repeats_pe_many_libs(k, gp, libes, scafolding_libes, output_dir, contigs_name, p);
}

//TODO: use only this one
void resolve_repeats_pe(size_t k, conj_graph_pack& gp, vector<PairedIndexT>& paired_index,
        vector<PairedIndexT>& mate_pair_index,
        vector<PairedIndexT>& scaffolder_index,
        vector<PairedIndexT>& mate_pair_scaffolding_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

	PairedIndexT paired_index_not_clust(gp.g);
	io::ReadStreamVector<SequencePairedReadStream> paired_streams = paired_binary_readers(true, cfg::get().ds.IS());
	FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index, gp.kmer_mapper, paired_index_not_clust, paired_streams, gp.k_value);
	PairedInfoLibraries libs;
	for (size_t i = 0; i < paired_index.size(); ++i) {
		libs.push_back(
				new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(),
						paired_index[i]));
	}
	PairedInfoLibraries libs_mp;
	for (size_t i = 0; i < mate_pair_index.size(); ++i) {
		libs_mp.push_back(
				new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(),
						mate_pair_index[i], true));
	}

    PairedInfoLibraries scaf_libs;
    for (size_t i = 0; i < scaffolder_index.size(); ++i){
    	scaf_libs.push_back(new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(),
    			scaffolder_index[i]));
    }
    PairedInfoLibraries scaf_libs_mp;
    for (size_t i = 0; i < mate_pair_scaffolding_index.size(); ++i){
    	scaf_libs_mp.push_back(new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(),
    			mate_pair_scaffolding_index[i], true));
    }

    vector<PairedInfoLibraries> libes;
    libes.push_back(libs);
    libes.push_back(libs_mp);
    vector<PairedInfoLibraries> scafolding_libes;
    scafolding_libes.push_back(scaf_libs);
    scafolding_libes.push_back(scaf_libs_mp);

    resolve_repeats_pe_many_libs(k, gp, libes, scafolding_libes,
          		output_dir, contigs_name, p);
}

void resolve_repeats_pe(size_t k, conj_graph_pack& gp, PairedInfoIndexT<Graph>& paired_index, PairedInfoIndexT<Graph>& scaffolder_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

	SingleThresholdFinder finder(cfg::get().ds.IS() - 2 * (cfg::get().ds.is_var()), cfg::get().ds.IS() + 2* (cfg::get().ds.is_var()), 100	);
	double threshold = finder.find_threshold();
	INFO("we found single threshold!! It is " << threshold);
	PairedIndexT paired_index_not_clust(gp.g);
	io::ReadStreamVector<SequencePairedReadStream> paired_streams = paired_binary_readers(true, cfg::get().ds.IS());
	FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index, gp.kmer_mapper, paired_index_not_clust, paired_streams, gp.k_value);
	PairedInfoLibraries libs;
	PairedInfoLibrary* lib = new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(), paired_index, paired_index_not_clust);
	lib->SetSingleThreshold(threshold);
	libs.push_back(lib);

    PairedInfoLibraries scaf_libs;
    scaf_libs.push_back(new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(), scaffolder_index));

    //resolve_repeats_pe_wcontigs(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
    vector<PairedInfoLibraries> libes;
    libes.push_back(libs);
    vector<PairedInfoLibraries> scafolding_libes;
    scafolding_libes.push_back(scaf_libs);

        //resolve_repeats_pe_wcontigs(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
     resolve_repeats_pe_many_libs(k, gp, libes, scafolding_libes,
          		output_dir, contigs_name, p);
}

/*void resolve_repeats_pe_mp(size_t k, conj_graph_pack& gp, PairedIndexT& paired_index,
        PairedIndexT& jpi,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(), paired_index));
    libs.push_back(new PairedInfoLibrary(k, gp.g, 93, 3280, 1500,  jpi, true));

    PairedInfoLibraries scaf_libs;


    resolve_repeats_pe(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
}*/

void resolve_repeats_pe_mp(size_t k, conj_graph_pack& gp, PairedIndexT& paired_index,
        PairedIndexT& mate_pair_index,
        PairedIndexT& scaffolder_index,
        PairedIndexT& mate_pair_scaffolding_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {
    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(), paired_index, scaffolder_index));

    PairedInfoLibraries mp_libs;
    mp_libs.push_back(new PairedInfoLibrary(k, gp.g, 100, 2200, 1000, mate_pair_index, mate_pair_scaffolding_index, true));

    PairedInfoLibraries scaf_libs;
    scaf_libs.push_back(new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(), scaffolder_index));

    PairedInfoLibraries mp_scaf_libs;
    mp_scaf_libs.push_back(new PairedInfoLibrary(k, gp.g, 93, 2200, 1000,  mate_pair_scaffolding_index, true));

    resolve_repeats_pe_mp(k, gp, libs, scaf_libs, mp_libs, mp_scaf_libs, output_dir, contigs_name, p);
}

void resolve_repeats_pe_mp(size_t k, conj_graph_pack& gp, PairedIndexT& paired_index,
        PairedIndexT& mate_pair_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, cfg::get().ds.RL(), cfg::get().ds.IS(), cfg::get().ds.is_var(), paired_index));

    PairedInfoLibraries mp_libs;
    mp_libs.push_back(new PairedInfoLibrary(k, gp.g, 100, 2200, 1000, mate_pair_index, true));

    PairedInfoLibraries scaf_libs;
    PairedInfoLibraries mp_scaf_libs;

    resolve_repeats_pe_mp(k, gp, libs, scaf_libs, mp_libs, mp_scaf_libs, output_dir, contigs_name, p);
}


} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
