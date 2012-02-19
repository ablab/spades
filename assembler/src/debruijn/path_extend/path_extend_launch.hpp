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
	params = p;
    if (!cfg::get().run_mode) {
        params.output.DisableAll();
        params.viz.DisableAll();
    }


    INFO("Path extend repeat resolving tool started");
    INFO("Using " << params.param_set.metric << " metric");
    INFO("Scaffolder is " << (cfg::get().pe_params.param_set.scaffolder_options.on ? "on" : "off"));

    PathInfoWriter path_writer;
	PathVisualizer visualizer(k);
	ContigWriter writer(gp.g, k);
	writer.writeEdges(etcDir + "before_resolve.fasta");

	INFO("Initializing weight counters");
	WeightCounter * wc = 0;
    WeightCounter * scaf_wc = 0;

    scaf_wc = new ReadCountWeightCounter(gp.g, scafolding_libs, params.param_set.extension_options.select_options.weight_threshold);
	if (params.param_set.metric == "read_count") {
	    wc = new ReadCountWeightCounter(gp.g, libs, params.param_set.extension_options.select_options.weight_threshold);
	}
	else if (params.param_set.metric == "path_cover") {
	    wc = new PathCoverWeightCounter(gp.g, libs, params.param_set.extension_options.select_options.weight_threshold, *params.param_set.extension_options.select_options.single_threshold);
	} else {
	    WARN("Unknown metric");
	    return;
	}

	wc->setNormalizeWeight(params.param_set.normalize_weight);
	wc->setNormalizeWightByCoverage(params.param_set.normalize_by_coverage);

	INFO("Initializing resolver");
	ExtensionChooser * seedEC;
	if (params.param_set.seed_selection.check_trusted) {
	    seedEC = new TrivialExtensionChooserWithPI(gp.g, wc);
	} else {
	    seedEC = new TrivialExtensionChooser(gp.g);
	}

	PathExtender * seedPE = new SimplePathExtender(gp.g, seedEC);
	seedPE->setMaxLoops(params.param_set.loop_removal.max_loops);

	PathExtendResolver resolver(gp.g, k);

	INFO("Finding seeds");
	auto seeds = resolver.makeSeeds(*seedPE, params.param_set.seed_selection.start_egde_coverage);
	INFO("Found " << seeds.size() << " seeds");

//	((SimplePathExtender*) seedPE)->getCoverageMap().printUncovered();
//	seeds.checkSymmetry();

	if (params.output.write_seeds) {
	    writer.writePaths(seeds, etcDir + "seeds.fasta");
	}
    if (params.viz.print_seeds) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "seeds.dot", "Seeds", seeds);
        path_writer.writePaths(seeds, etcDir + "seeds.data");
    }

    if ( params.param_set.seed_selection.min_coverage > 0) {
        INFO("Filtering seeds by coverage " << params.param_set.seed_selection.min_coverage);
        CoveragePathFilter coverageFilter(gp.g, params.param_set.seed_selection.min_coverage);
        coverageFilter.filter(seeds);
    }

    SimpleExtensionChooser * extensionChooser = new SimpleExtensionChooser(gp.g, wc, params.param_set.extension_options.select_options.priority_coeff);
    ScaffoldingExtensionChooser * scafExtensionChooser = new ScaffoldingExtensionChooser(gp.g, scaf_wc, params.param_set.extension_options.select_options.priority_coeff);
	ScaffoldingPathExtender * mainPE = new ScaffoldingPathExtender(gp.g, extensionChooser, scafExtensionChooser);

	seeds.SortByLength();
	INFO("Extending seeds");

	auto paths = resolver.extendSeeds(seeds, *mainPE);
    if (params.output.write_overlaped_paths) {
        writer.writePaths(paths, etcDir + "overlaped_paths.fasta");
    }
    if (params.viz.print_overlaped_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "overlaped_paths.dot", "Overlaped_paths", paths);
        path_writer.writePaths(paths, etcDir + "overlaped_paths.data");
    }

//    mainPE->getCoverageMap().printUncovered();
//    paths.checkSymmetry();

	if (params.param_set.filter_options.remove_overlaps) {
	    INFO("Removing overlaps");
	    resolver.removeOverlaps(paths, mainPE->GetCoverageMap());
	}

	mainPE->GetCoverageMap().PrintUncovered();
	paths.CheckSymmetry();

	INFO("Adding uncovered edges");
	resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());

	paths.SortByLength();

//	mainPE->getCoverageMap().printUncovered();
//	paths.checkSymmetry();

	INFO("Start validating");
    PathValidator pathValidator(gp.g, wc, scaf_wc);
    vector<bool> validationResult = pathValidator.ValidatePaths(paths);

	INFO("Found " << paths.size() << " contigs");
	writer.writePaths(paths, output_dir + contigs_name);
    if (params.output.write_paths) {
        writer.writePathEdges(paths, output_dir + "paths.dat");
    }

    if (params.viz.print_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "paths.dot", "Final_paths", paths);
        path_writer.writePaths(paths, etcDir + "paths.data");
    }
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


void resolve_repeats_pe(size_t k, conj_graph_pack& gp, PairedInfoIndex<Graph>& paired_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, paired_index));

    PairedInfoLibraries scaf_libs;

    resolve_repeats_pe(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
}



void resolve_repeats_pe(size_t k, conj_graph_pack& gp, PairedInfoIndex<Graph>& paired_index, PairedInfoIndex<Graph>& scaffolder_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {


    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, paired_index));

    PairedInfoLibraries scaf_libs;
    scaf_libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, scaffolder_index));

    resolve_repeats_pe(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
}

void resolve_repeats_pe_wj(size_t k, conj_graph_pack& gp, PairedInfoIndex<Graph>& paired_index,
        PairedInfoIndex<Graph>& jpi,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;
    libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, paired_index));
    libs.push_back(new PairedInfoLibrary(k, gp.g, 150, 7495, 1290, jpi));

    PairedInfoLibraries scaf_libs;

    resolve_repeats_pe(k, gp, libs, scaf_libs, output_dir, contigs_name, p);
}

void resolve_repeats_pe_wj(size_t k, conj_graph_pack& gp, PairedInfoIndex<Graph>& paired_index,
        PairedInfoIndex<Graph>& jpi,
        PairedInfoIndex<Graph>& scaffolder_index,
        PairedInfoIndex<Graph>& jspi,
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
    if (params.debug_output) {
        make_dir(etcDir);
    }

    params = p;
    if (!cfg::get().run_mode) {
        params.output.DisableAll();
        params.viz.DisableAll();
    }

    INFO("Path extend scaffolding tool started");

    PathInfoWriter path_writer;
    PathVisualizer visualizer(k);
    ContigWriter writer(gp.g, k);

    if (params.debug_output) {
        writer.writeEdges(etcDir + "before_resolve.fasta");
    }

    WeightCounter * scaf_wc = 0;
    scaf_wc = new ReadCountWeightCounter(gp.g, scafolding_libs, params.param_set.extension_options.select_options.weight_threshold);
    scaf_wc->setNormalizeWeight(false);
    scaf_wc->setNormalizeWightByCoverage(false);

    PathExtendResolver resolver(gp.g, k);

    INFO("Finding seeds");
    auto seeds = resolver.ReturnEdges(params.param_set.seed_selection.start_egde_coverage);
    INFO("Found " << seeds.size() << " seeds");

//  ((SimplePathExtender*) seedPE)->getCoverageMap().printUncovered();
//  seeds.checkSymmetry();

    if (params.output.write_seeds) {
        writer.writePaths(seeds, etcDir + "seeds.fasta");
    }
    if (params.viz.print_seeds) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "seeds.dot", "Seeds", seeds);
        path_writer.writePaths(seeds, etcDir + "seeds.data");
    }

    if ( params.param_set.seed_selection.min_coverage > 0) {
        INFO("Filtering seeds by coverage " << params.param_set.seed_selection.min_coverage);
        CoveragePathFilter coverageFilter(gp.g, params.param_set.seed_selection.min_coverage);
        coverageFilter.filter(seeds);
    }


    ScaffoldingExtensionChooser * scafExtensionChooser = new ScaffoldingExtensionChooser(gp.g, scaf_wc, params.param_set.extension_options.select_options.priority_coeff);

    ScaffoldingOnlyPathExtender * mainPE = new ScaffoldingOnlyPathExtender(gp.g, scafExtensionChooser);

    mainPE->setMaxLoops(params.param_set.loop_removal.max_loops);
    mainPE->setInvestigateShortLoops(params.param_set.loop_removal.inspect_short_loops);
    seeds.SortByLength();
    INFO("Extending seeds");
    auto paths = resolver.extendSeeds(seeds, *mainPE);


    if (params.output.write_overlaped_paths) {
        writer.writePaths(paths, etcDir + "overlaped_paths.fasta");
    }
    if (params.viz.print_overlaped_paths) {
        visualizer.writeGraphWithPathsSimple(gp, etcDir + "overlaped_paths.dot", "Overlaped_paths", paths);
        path_writer.writePaths(paths, etcDir + "overlaped_paths.data");
    }

//    mainPE->getCoverageMap().printUncovered();
//    paths.checkSymmetry();

    if (params.param_set.filter_options.remove_overlaps) {
        INFO("Removing overlaps");
        resolver.removeOverlaps(paths, mainPE->GetCoverageMap());
    }

    mainPE->GetCoverageMap().PrintUncovered();
    paths.CheckSymmetry();

    INFO("Adding uncovered edges");
    resolver.addUncoveredEdges(paths, mainPE->GetCoverageMap());

    paths.SortByLength();

    INFO("Found " << paths.size() << " scaffolds");
    writer.writePaths(paths, output_dir + contigs_name);
    if (params.output.write_paths) {
        writer.writePathEdges(paths, output_dir + "paths.dat");
    }

    INFO("Path extend scaffolding tool finished");
}



void scaffold_pe(size_t k, conj_graph_pack& gp, PairedInfoIndex<Graph>& scaffolder_index,
        const std::string& output_dir, const std::string& contigs_name, const pe_config::MainPEParamsT& p) {

    PairedInfoLibraries libs;

    PairedInfoLibraries scaf_libs;
    scaf_libs.push_back(new PairedInfoLibrary(k, gp.g, *cfg::get().ds.RL, *cfg::get().ds.IS, *cfg::get().ds.is_var, scaffolder_index));

    scaffold_pe(k, gp, scaf_libs, output_dir, contigs_name, p);
}

} /* path_extend */



#endif /* PATH_EXTEND_LAUNCH_HPP_ */
