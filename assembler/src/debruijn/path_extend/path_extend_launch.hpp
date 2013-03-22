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

//#include "gap_jumper.hpp"

#include "../include/omni/edges_position_handler.hpp"

namespace path_extend {

using namespace debruijn_graph;


void resolve_repeats_pe(size_t k, conj_graph_pack& gp, PairedInfoLibraries& libs, PairedInfoLibraries& scafolding_libs,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {
	make_dir(output_dir);
	std::string etcDir = output_dir + "path_extend/";
	make_dir(etcDir);
  if (!cfg::get().run_mode) {
    cfg::get_writable().pe_params.output.DisableAll();
    cfg::get_writable().pe_params.viz.DisableAll();
  }

  INFO("Path extend repeat resolving tool started");
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
	    seedEC = new TrivialExtensionChooserWithPI(gp.g, wc);
	} else {
	    seedEC = new TrivialExtensionChooser(gp.g);
	}

	PathExtender * seedPE = new SimplePathExtender(gp.g, seedEC);
	seedPE->setMaxLoops(cfg::get().pe_params.param_set.loop_removal.max_loops);

	PathExtendResolver resolver(gp.g, k);

	INFO("Finding seeds");
	auto seeds = resolver.makeSeeds(*seedPE, cfg::get().pe_params.param_set.seed_selection.start_egde_coverage);
	INFO("Found " << seeds.size() << " seeds");

//	((SimplePathExtender*) seedPE)->getCoverageMap().printUncovered();
//	seeds.checkSymmetry();

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
    ScaffoldingExtensionChooser * scafExtensionChooser = new ScaffoldingExtensionChooser(gp.g, scaf_wc, cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);
	ScaffoldingPathExtender * mainPE = new ScaffoldingPathExtender(gp.g, extensionChooser, scafExtensionChooser);

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

//    mainPE->getCoverageMap().printUncovered();
//    paths.checkSymmetry();

	if (cfg::get().pe_params.param_set.filter_options.remove_overlaps) {
	    INFO("Removing overlaps");
	    resolver.removeOverlaps(paths, mainPE->GetCoverageMap());
	}

	//mainPE->GetCoverageMap().PrintUncovered();
	//paths.CheckSymmetry();

	DEBUG("Adding uncovered edges");
	resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());

	paths.SortByLength();

//	mainPE->getCoverageMap().printUncovered();
//	paths.checkSymmetry();

//	DEBUG("Start validating");
//    PathValidator pathValidator(gp.g, wc, scaf_wc);
//    vector<bool> validationResult = pathValidator.ValidatePaths(paths);

	INFO("Found " << paths.size() << " contigs");
	writer.writePaths(paths, output_dir + contigs_name);
    if (cfg::get().pe_params.output.write_paths) {
        writer.writePathEdges(paths, output_dir + "paths.dat");
    }

    if (cfg::get().pe_params.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "paths.dot", "Final_paths", paths);
        path_writer.writePaths(paths, etcDir + "paths.data");
    }
    writer.writePaths(paths, output_dir + contigs_name);

//    PathContainer brokenScaffolds;
//    for (size_t i = 0; i < paths.size(); ++i) {
//        BidirectionalPath * path = paths.Get(i);
//        size_t last_pos = 0;
//        for (size_t j = 0; j < path->Size(); ++j) {
//            if (path->GapAt(j) > (int) gp.g.k()) {
//                BidirectionalPath * p = new BidirectionalPath(path->SubPath(last_pos, j));
//                BidirectionalPath * cp = new BidirectionalPath(p->Conjugate());
//                last_pos = j;
//                brokenScaffolds.AddPair(p, cp);
//            }
//        }
//        BidirectionalPath * p = new BidirectionalPath(path->SubPath(last_pos));
//        BidirectionalPath * cp = new BidirectionalPath(p->Conjugate());
//        brokenScaffolds.AddPair(p, cp);
//    }
//    writer.writePaths(brokenScaffolds, output_dir + "broken_scaffolds.fasta");

//    size_t count = 0;
//    ofstream out("./sinks.log");
//    for (size_t i = 0; i < paths.size(); ++ i) {
//    	 if (gp.edge_pos.IsConsistentWithGenome(paths.Get(i)->ToVector())) {
//    		 if (gp.g.OutgoingEdgeCount(gp.g.EdgeEnd(paths.Get(i)->Back())) == 0) {
//                out << paths.Get(i)->GetId() << " " << gp.g.int_id(paths.Get(i)->Back()) << "\n";
//    		 }
//    		 count ++;
//    	 }
//    }
//    INFO("Number of paths consistent with genome: " << count);
//
//    INFO("Total number of paths: " << paths.size());

	INFO("Path extend repeat resolving tool finished");

}


void resolve_repeats_pe(size_t k, conj_graph_pack& gp, PairedInfoIndexT<Graph>& paired_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, paired_index));

    PairedInfoLibraries scaf_libs;

    resolve_repeats_pe(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
}



void resolve_repeats_pe(size_t k, conj_graph_pack& gp, PairedInfoIndexT<Graph>& paired_index, PairedInfoIndexT<Graph>& scaffolder_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {


    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, paired_index));

    PairedInfoLibraries scaf_libs;
    scaf_libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, scaffolder_index));

    resolve_repeats_pe(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
}

void resolve_repeats_pe_wj(size_t k, conj_graph_pack& gp, PairedInfoIndexT<Graph>& paired_index,
        PairedInfoIndexT<Graph>& jpi,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, paired_index));
    libs.push_back(new PairedInfoLibrary(k, gp.g, 150, 7495, 1290, jpi));

    PairedInfoLibraries scaf_libs;

    resolve_repeats_pe(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
}

void resolve_repeats_pe_wj(size_t k, conj_graph_pack& gp, PairedInfoIndexT<Graph>& paired_index,
        PairedInfoIndexT<Graph>& jpi,
        PairedInfoIndexT<Graph>& scaffolder_index,
        PairedInfoIndexT<Graph>& jspi,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, paired_index));
    libs.push_back(new PairedInfoLibrary(k, gp.g, 150, 7495, 1290, jpi));

    PairedInfoLibraries scaf_libs;
    scaf_libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, scaffolder_index));
    libs.push_back(new PairedInfoLibrary(k, gp.g, 150, 7495, 1290, jspi));

    resolve_repeats_pe(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
}




void scaffold_pe(size_t k, conj_graph_pack& gp, PairedInfoLibraries& scafolding_libs,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {
    make_dir(output_dir);

    std::string etcDir = output_dir + "path_extend/";
    if (cfg::get().pe_params.debug_output) {
        make_dir(etcDir);
    }

    if (!cfg::get().run_mode) {
        cfg::get_writable().pe_params.output.DisableAll();
        cfg::get_writable().pe_params.viz.DisableAll();
    }

    INFO("Path extend scaffolding tool started");

    PathInfoWriter path_writer;
    PathVisualizer visualizer(k);
    ContigWriter writer(gp, k);

    if (cfg::get().pe_params.debug_output) {
        writer.writeEdges(etcDir + "before_resolve.fasta");
        visualizer.writeGraphSimple(gp, etcDir + "before_resolve.dot", "before_resolve");
    }

    WeightCounter * scaf_wc = 0;
    scaf_wc = new ReadCountWeightCounter(gp.g, scafolding_libs, cfg::get().pe_params.param_set.extension_options.select_options.weight_threshold);
    scaf_wc->setNormalizeWeight(false);
    scaf_wc->setNormalizeWightByCoverage(false);

    PathExtendResolver resolver(gp.g, k);

    INFO("Finding seeds");
    auto seeds = resolver.ReturnEdges(cfg::get().pe_params.param_set.seed_selection.start_egde_coverage);
    INFO("Found " << seeds.size() << " seeds");

//  ((SimplePathExtender*) seedPE)->getCoverageMap().printUncovered();
//  seeds.checkSymmetry();

    if (cfg::get().pe_params.output.write_seeds) {
        writer.writePaths(seeds, etcDir + "seeds.fasta");
    }
    if (cfg::get().pe_params.viz.print_seeds) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "seeds.dot", "Seeds", seeds);
        path_writer.writePaths(seeds, etcDir + "seeds.data");
    }

    if ( cfg::get().pe_params.param_set.seed_selection.min_coverage > 0) {
        INFO("Filtering seeds by coverage " << cfg::get().pe_params.param_set.seed_selection.min_coverage);
        CoveragePathFilter coverageFilter(gp.g, cfg::get().pe_params.param_set.seed_selection.min_coverage);
        coverageFilter.filter(seeds);
    }


    ScaffoldingExtensionChooser * scafExtensionChooser = new ScaffoldingExtensionChooser(gp.g, scaf_wc, cfg::get().pe_params.param_set.extension_options.select_options.priority_coeff);

    ScaffoldingOnlyPathExtender * mainPE = new ScaffoldingOnlyPathExtender(gp.g, scafExtensionChooser);

    mainPE->setMaxLoops(cfg::get().pe_params.param_set.loop_removal.max_loops);
    mainPE->setInvestigateShortLoops(cfg::get().pe_params.param_set.loop_removal.inspect_short_loops);
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

//    mainPE->getCoverageMap().printUncovered();
//    paths.checkSymmetry();

    if (cfg::get().pe_params.param_set.filter_options.remove_overlaps) {
        INFO("Removing overlaps");
        resolver.removeOverlaps(paths, mainPE->GetCoverageMap());
    }

    mainPE->GetCoverageMap().PrintUncovered();
    paths.CheckSymmetry();

    DEBUG("Adding uncovered edges");
    resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());

    paths.SortByLength();

    INFO("Found " << paths.size() << " scaffolds");
    writer.writePaths(paths, output_dir + contigs_name);
    if (cfg::get().pe_params.output.write_paths) {
        writer.writePathEdges(paths, output_dir + "paths.dat");
    }

    INFO("Path extend scaffolding tool finished");
}



void scaffold_pe(size_t k, conj_graph_pack& gp, PairedInfoIndexT<Graph>& scaffolder_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;

    PairedInfoLibraries scaf_libs;
    scaf_libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, scaffolder_index));

    scaffold_pe(k, gp, scaf_libs, output_dir, contigs_name, p);
}

} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
