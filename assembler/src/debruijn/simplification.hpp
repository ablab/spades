/*
 * simplification.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"
#include "construction.hpp"
#include "gap_closer.hpp"
#include "omni_labelers.hpp"
#include "omni/omni_tools.hpp"

namespace debruijn_graph {
void simplify_graph(PairedReadStream& stream, conj_graph_pack& gp,
		paired_info_index& paired_index);
} // debruijn_graph

// move impl to *.cpp
namespace debruijn_graph {

template<size_t k>
void PrintWeightDistribution(Graph &g, const string &file_name) {
	ofstream os(file_name.c_str());
	for(auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		vector<EdgeId> v1 = g.OutgoingEdges(g.EdgeStart(*it));
		vector<EdgeId> v2 = g.IncomingEdges(g.EdgeEnd(*it));
		bool eq = false;
		if(v1.size() == 2 && v2.size() == 2)
			if((v1[0] == v2[0] && v1[1] == v2[1]) || (v1[0] == v2[1] && v1[0] == v2[1]))
				eq = false;
		if(g.length(*it) > k - 10 && g.length(*it) <= k + 1 && g.OutgoingEdgeCount(g.EdgeStart(*it))>= 2 && g.IncomingEdgeCount(g.EdgeEnd(*it))>= 2 && !eq)
			os << g.coverage(*it) << endl;
	}
	os.close();
}

void simplify_graph(conj_graph_pack& gp) {
	using namespace omnigraph;

	exec_construction(gp);

	INFO("STAGE == Simplifying graph");

//	PrintWeightDistribution<K>(gp.g, "distribution.txt");
	EdgeQuality<Graph> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
	total_labeler_graph_struct graph_struct(gp.g, &gp.int_ids, &gp.edge_pos);
	total_labeler tot_lab(&graph_struct);

	CompositeLabeler<Graph> labeler(tot_lab, edge_qual);

	detail_info_printer printer(gp, labeler, cfg::get().output_dir, "graph.dot");
    printer(ipp_before_first_gap_closer);
	
    if (cfg::get().gap_closer_enable && cfg::get().gc.before_simplify)
		CloseGap<K>(gp, cfg::get().gc.use_extended_mapper);
	

//	QualityLoggingRemovalHandler<Graph> qual_removal_handler(gp.g, edge_qual);
	QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g,
			edge_qual,
			labeler,
			cfg::get().output_dir);

	boost::function<void(EdgeId)> removal_handler_f = boost::bind(
//			&QualityLoggingRemovalHandler<Graph>::HandleDelete,
			&QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
			boost::ref(qual_removal_handler), _1);


	SimplifyGraph(gp, removal_handler_f, labeler, printer, 10
			/*, etalon_paired_index*/);

	AvgCovereageCounter<Graph> cov_counter(gp.g);
	cfg::get_writable().ds.avg_coverage = cov_counter.Count();

	if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
		CloseGap<K>(gp, cfg::get().gc.use_extended_mapper);
	
	//  ProduceInfo<k>(g, index, *totLab, genome, output_folder + "simplified_graph.dot", "simplified_graph");


	//experimental
//	if (cfg::get().paired_mode) {
//		INFO("Pair info aware ErroneousConnectionsRemoval");
//		RemoveEroneousEdgesUsingPairedInfo(gp.g, paired_index);
//		INFO("Pair info aware ErroneousConnectionsRemoval stats");
//		CountStats<K>(gp.g, gp.index, gp.genome);
//	}
	//experimental

	//	ProduceDetailedInfo<k>(g, index, labeler, genome, output_folder + "with_pair_info_edges_removed/",	"graph.dot", "no_erroneous_edges_graph");

	//  WriteGraphComponents<k>(g, index, *totLab, genome, output_folder + "graph_components" + "/", "graph.dot",
	//            "graph_component", cfg::get().ds.IS);

	//  number_of_components = PrintGraphComponents(output_folder + "graph_components/graph", g,
	//            cfg::get().ds.IS, int_ids, paired_index, EdgePos);
}

void load_simplification(conj_graph_pack& gp, files_t* used_files) {
	fs::path p = fs::path(cfg::get().load_from) / "simplified_graph";
	used_files->push_back(p);

	ScanGraphPack(p.string(), gp);
	load_estimated_params(p.string());
}

void save_simplification(conj_graph_pack& gp) {
	if (cfg::get().make_saves) {
		fs::path p = fs::path(cfg::get().output_saves) / "simplified_graph";
		PrintGraphPack(p.string(), gp);
		write_estimated_params(p.string());
	}

	//todo temporary solution!!!
	OutputContigs(gp.g, cfg::get().additional_contigs);
	OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");
    cfg::get_writable().final_contigs_file = cfg::get().output_dir + "final_contigs.fasta";

// run script automatically takes simplified contigs from correct path

//	OutputContigs(gp.g,
//			cfg::get().output_root + "../" + cfg::get().additional_contigs);
}

void exec_simplification(conj_graph_pack& gp) {
	if (cfg::get().entry_point <= ws_simplification) {
		simplify_graph(gp);
		save_simplification(gp);
	} else {
		INFO("Loading Simplification");

		files_t used_files;
		load_simplification(gp, &used_files);
		link_files_by_prefix(used_files, cfg::get().output_saves);
	}
}

} //debruijn_graph
