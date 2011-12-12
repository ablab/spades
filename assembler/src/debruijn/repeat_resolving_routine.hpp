/*
 * repeat_resolving_routine.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"

#include "logging.hpp"
#include "repeat_resolving.hpp"
#include "distance_estimation_routine.hpp"
#include "path_set_graph_constructor.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"
#include "resolved_pair_info.hpp"
#include "debruijn_stats.hpp"
//typedef io::IReader<io::SingleRead> ReadStream;
//typedef io::IReader<io::PairedRead> PairedReadStream;
////typedef io::RCReaderWrapper<io::SingleRead> RCStream;
//typedef io::MultifileReader<io::SingleRead> MultiFileStream;
typedef io::CarefulFilteringReaderWrapper<io::SingleRead> CarefulFilteringStream;

namespace debruijn_graph

{
void resolve_repeats(PairedReadStream& stream, const Sequence& genome);
} // debruijn_graph

// move impl to *.cpp

namespace debruijn_graph {

void FillContigNumbers(map<ConjugateDeBruijnGraph::EdgeId, int>& contigNumbers
		, ConjugateDeBruijnGraph& cur_graph) {
	int cur_num = 0;
	set<ConjugateDeBruijnGraph::EdgeId> edges;
	for (auto iter = cur_graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		if (edges.find(*iter) == edges.end()) {
			contigNumbers[*iter] = cur_num;
			cur_num++;
			edges.insert(cur_graph.conjugate(*iter));
		}
	}
}

void FillContigNumbers(
		map<NonconjugateDeBruijnGraph::EdgeId, int>& contigNumbers
		, NonconjugateDeBruijnGraph& cur_graph) {
	int cur_num = 0;
	for (auto iter = cur_graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		contigNumbers[*iter] = cur_num;
		cur_num++;
	}
}

int ContigNumber(map<NonconjugateDeBruijnGraph::EdgeId, int>& contigNumbers , NonconjugateDeBruijnGraph::EdgeId eid,  NonconjugateDeBruijnGraph& cur_graph){
	if (contigNumbers.find(eid) != contigNumbers.end())
		return(contigNumbers[eid]);
	else {
		WARN("Deleted edge");
		return -1;
	}
}

int ContigNumber(map<ConjugateDeBruijnGraph::EdgeId, int>& contigNumbers , ConjugateDeBruijnGraph::EdgeId eid,  ConjugateDeBruijnGraph& cur_graph){
	if (contigNumbers.find(eid) != contigNumbers.end())
		return(contigNumbers[eid]);
	else
		if (contigNumbers.find(cur_graph.conjugate(eid)) != contigNumbers.end())
			return(contigNumbers[cur_graph.conjugate(eid)]);
	else {
		WARN("Deleted edge");
		return -1;
	}
}

template<size_t k, class graph_pack>
void SelectReadsForConsensusBefore(graph_pack& etalon_gp, typename graph_pack::graph_t& cur_graph,
        EdgeLabelHandler<typename graph_pack::graph_t>& LabelsAfter,
        const EdgeIndex<K + 1, typename graph_pack::graph_t>& index ,vector<ReadStream *>& reads
        , string& consensus_output_dir)
{
    INFO("ReadMapping started");
    map<typename graph_pack::graph_t::EdgeId, int> contigNumbers;
    int cur_num = 0;
    FillContigNumbers(contigNumbers, cur_graph);
    for(auto iter = etalon_gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter){
    	DEBUG("Edge number:" << etalon_gp.int_ids.ReturnIntId(*iter) << " is contained in contigs" );
        set<typename graph_pack::graph_t::EdgeId> images =
                            LabelsAfter.edge_inclusions[*iter];
        for (auto it = images.begin(); it != images.end(); ++it) {
        	DEBUG(ContigNumber(contigNumbers, *it, cur_graph) << ", ");
        }
    }
    cur_num = contigNumbers.size();
    INFO(cur_num << "contigs");
    for (int i = 1; i < 3; i++) {
        int read_num = 0;
        osequencestream* mapped_reads[5000];
        for (int j = 0; j < cur_num; j++) {
            string output_filename = consensus_output_dir + ToString(j)
                    + "_reads" + ToString(i) + ".fa";
            osequencestream* tmp = new osequencestream(output_filename);
//          mapped_reads.push_back(tmp);
            mapped_reads[j] = tmp;
        }
        SingleReadMapper<k, typename graph_pack::graph_t> rm(etalon_gp.g, index);
        INFO("mapping reads from pair"<< i);
        while (!reads[i - 1]->eof()) {
            io::SingleRead cur_read;

            (* reads[i - 1]) >> cur_read;
            vector<typename graph_pack::graph_t::EdgeId> res = rm.GetContainingEdges(
                    cur_read);
            read_num++;
            TRACE(
                    read_num<< " mapped to"<< res.size() <<" contigs :, read"<< cur_read.sequence());
//          map_quantity += res.size();
            for (size_t ii = 0; ii < res.size(); ii++) {
                TRACE("counting number "<< contigNumbers[res[ii]]);
                	if (ContigNumber(contigNumbers, res[ii], cur_graph) != -1)
                		(*mapped_reads[ContigNumber(contigNumbers, res[ii], cur_graph)])
                			<< cur_read.sequence();
                	else
                		WARN("No edges containing" <<etalon_gp.int_ids.ReturnIntId(res[ii]));
            }
        }
    }
}


template<size_t k, class graph_pack>
void SelectReadsForConsensus(graph_pack& etalon_gp, typename graph_pack::graph_t& cur_graph,
        EdgeLabelHandler<typename graph_pack::graph_t>& LabelsAfter,
        const EdgeIndex<K + 1, typename graph_pack::graph_t>& index ,vector<ReadStream *>& reads
        , string& consensus_output_dir)
{
    INFO("ReadMapping started");
    map<typename graph_pack::graph_t::EdgeId, int> contigNumbers;
    int cur_num = 0;
    FillContigNumbers(contigNumbers, cur_graph);
    for(auto iter = etalon_gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter){
    	DEBUG("Edge number:" << etalon_gp.int_ids.ReturnIntId(*iter) << " is contained in contigs" );
        set<typename graph_pack::graph_t::EdgeId> images =
                            LabelsAfter.edge_inclusions[*iter];
        for (auto it = images.begin(); it != images.end(); ++it) {
        	DEBUG(ContigNumber(contigNumbers, *it, cur_graph) << ", ");
        }
    }
    cur_num = contigNumbers.size();
    INFO(cur_num << "contigs");
    for (int i = 1; i < 3; i++) {
        int read_num = 0;
        osequencestream* mapped_reads[5000];
        for (int j = 0; j < cur_num; j++) {
            string output_filename = consensus_output_dir + ToString(j)
                    + "_reads" + ToString(i) + ".fa";
            osequencestream* tmp = new osequencestream(output_filename);
//          mapped_reads.push_back(tmp);
            mapped_reads[j] = tmp;
        }
        SingleReadMapper<k, typename graph_pack::graph_t> rm(etalon_gp.g, index);
        INFO("mapping reads from pair"<< i);
        while (!reads[i - 1]->eof()) {
            io::SingleRead cur_read;

            (* reads[i - 1]) >> cur_read;
            vector<typename graph_pack::graph_t::EdgeId> res = rm.GetContainingEdges(
                    cur_read);
            read_num++;
            TRACE(
                    read_num<< " mapped to"<< res.size() <<" contigs :, read"<< cur_read.sequence());
//          map_quantity += res.size();
            for (size_t ii = 0; ii < res.size(); ii++) {
                TRACE("counting number "<< contigNumbers[res[ii]]);
                set<typename graph_pack::graph_t::EdgeId> images =
                        LabelsAfter.edge_inclusions[res[ii]];
                for (auto iter = images.begin(); iter != images.end(); ++iter)
                	if (ContigNumber(contigNumbers, *iter, cur_graph) != -1)
                		(*mapped_reads[ContigNumber(contigNumbers, *iter, cur_graph)])
                			<< cur_read.sequence();
                	else
                		WARN("No edges containing" <<etalon_gp.int_ids.ReturnIntId(res[ii]));
            }
        }
    }
}

template<class graph_pack>
void CleanIsolated(graph_pack& gp){
    for(auto iter = gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
    	typename graph_pack::graph_t::VertexId start, end;
    	start = gp.g.EdgeStart(*iter);
       	end = gp.g.EdgeEnd(*iter);
       	TRACE (gp.g.CheckUniqueOutgoingEdge(start)<<" "<<  gp.g.IsDeadStart(start) <<" "<< gp.g.CheckUniqueIncomingEdge(end) <<" "<<gp.g.IsDeadEnd(end));
       	if (gp.g.CheckUniqueOutgoingEdge(start) && gp.g.IsDeadStart(start) && gp.g.CheckUniqueIncomingEdge(end) && gp.g.IsDeadEnd(end) )
       		gp.g.DeleteEdge(*iter);
    }
}

string GeneratePostfix(){
	string s = "_";
	if (cfg::get().rr.symmetric_resolve)
		s += "sym_";
	else
		s += "nonsym_";

	if (cfg::get().advanced_estimator_mode)
		s += "advanced_est_";
	else
		s += "usual_est_";

	if (cfg::get().late_paired_info)
		s += "late_pi_";
	else
		s += "early_pi_";
	s += "k";
	s += ToString(K);
	if (cfg::get().path_set_graph){
		s += "_path_set";
	}
	if (cfg::get().rr.mode == 2) {
		s += "_mode2";
	} else {
		s += "_nv";
		s += ToString(cfg::get().rr.near_vertex);
	}
	s += ".fasta";
	return s;
}

template<class graph_pack>
void ProduceResolvedPairedInfo(graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp,
		EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
		PairedInfoIndex<typename graph_pack::graph_t>& resolved_graph_paired_info) {

    INFO("Generating paired info for resolved graph");
    ResolvedGraphPairInfoCounter<typename graph_pack::graph_t> resolved_graph_paired_info_counter(
    			origin_gp.g,
    			clustered_index,
    			resolved_gp.g,
    			labels_after);

    resolved_graph_paired_info_counter.FillResolvedGraphPairedInfo(resolved_graph_paired_info);
    INFO("Paired info for resolved graph generated");
}

template<class graph_pack>
void SaveResolvedPairedInfo(graph_pack& resolved_gp,
		PairedInfoIndex<typename graph_pack::graph_t> resolved_graph_paired_info,
		const string& graph_name, const string& subfolder) {

	std::string rr_filename = cfg::get().output_dir + subfolder + graph_name;
	INFO("Saving graph and paired info to " << rr_filename);
	PrintWithClusteredIndex(rr_filename, resolved_gp, resolved_graph_paired_info);
	INFO("Saved");
}

template<class graph_pack>
void process_resolve_repeats(graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp, const string& graph_name,
		const string& subfolder = "", bool output_contigs = true) {
	EdgeLabelHandler<typename graph_pack::graph_t> labels_after(resolved_gp.g,
			origin_gp.g);
	DEBUG("New index size: "<< clustered_index.size());
	// todo: make printGraph const to its arguments

	// todo: possibly we don't need it
//    if (cfg::get().rectangle_mode)
//        RectangleResolve(clustered_index, origin_gp.g, cfg::get().output_root + "tmp/", cfg::get().output_dir);
	string postfix = GeneratePostfix();
    typedef TotalLabelerGraphStruct<typename graph_pack::graph_t> total_labeler_gs;
    typedef TotalLabeler           <typename graph_pack::graph_t> total_labeler;
    total_labeler_gs graph_struct_before(origin_gp  .g, &origin_gp  .int_ids, &origin_gp  .edge_pos, NULL);
    total_labeler tot_labeler_before(&graph_struct_before);


	if (cfg::get().path_set_graph ) {
		INFO("testing path-set graphs");
		PathSetGraphConstructor<graph_pack> path_set_constructor(origin_gp, clustered_index,  resolved_gp);
		path_set_constructor.ConstructOverlap();
		unordered_map<typename graph_pack::graph_t::EdgeId, typename graph_pack::graph_t::EdgeId> edge_labels =
		path_set_constructor.GetEdgeLabels();
		labels_after.FillLabels(edge_labels);
		INFO("testing ended");
	}	else {

//    CleanIsolated(origin_gp);
		ResolveRepeats(origin_gp  .g, origin_gp  .int_ids, clustered_index, origin_gp  .edge_pos,
                   resolved_gp.g, resolved_gp.int_ids,                  resolved_gp.edge_pos,
                   cfg::get().output_dir + subfolder +"resolve_" + graph_name +  "/", labels_after);

    //Generating paired info for resolved graph
		PairedInfoIndex<typename graph_pack::graph_t> resolved_graph_paired_info(resolved_gp.g);
		ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp, labels_after, resolved_graph_paired_info);
		SaveResolvedPairedInfo(resolved_gp, resolved_graph_paired_info, graph_name + "_resolved", subfolder);
    //Paired info for resolved graph generated
	}
    if (output_contigs) {

       	OutputContigs(resolved_gp.g, cfg::get().output_dir + "after_rr_before_simplify" + postfix);
    	OutputContigs(origin_gp.g, cfg::get().output_dir + "before_resolve" + postfix);
    }
    INFO("Total labeler start");
    total_labeler_gs graph_struct_after (resolved_gp.g, &resolved_gp.int_ids, &resolved_gp.edge_pos, &labels_after);
    total_labeler tot_labeler_after(&graph_struct_after, &graph_struct_before);
    omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after, cfg::get().output_dir + subfolder + graph_name + "_3_resolved.dot", "no_repeat_graph");

    INFO("Total labeler finished");

    INFO("---Clearing resolved graph---");



    EdgeRemover<typename graph_pack::graph_t> edge_remover(resolved_gp.g, false);
    for (int i = 0; i < 3; ++i)
    {
        ClipTipsForResolver(resolved_gp.g);
//        BulgeRemoveWrap      (resolved_gp.g);
        RemoveLowCoverageEdges(resolved_gp.g, edge_remover, 3, i);

//        RemoveRelativelyLowCoverageEdges(resolved_gp.g);
        omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,

        	cfg::get().output_dir + subfolder + ToString(i) +  "b_4_cleared.dot",
        			"no_repeat_graph");
	}

	INFO("---Cleared---");

    //Generating paired info for resolved graph
    PairedInfoIndex<typename graph_pack::graph_t> resolved_cleared_graph_paired_info(resolved_gp.g);
    if (cfg::get().path_set_graph == false) {
    	ProduceResolvedPairedInfo(origin_gp, clustered_index, resolved_gp, labels_after, resolved_cleared_graph_paired_info);
    	SaveResolvedPairedInfo(resolved_gp, resolved_cleared_graph_paired_info, graph_name + "_resolved_cleared", subfolder);
	}
    //Paired info for resolved graph generated

	INFO("---Output Contigs---");

	if (output_contigs)
		OutputContigs(resolved_gp.g,
				cfg::get().output_dir + "resolved_and_cleared" + postfix);

	omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,

	cfg::get().output_dir + subfolder + graph_name + "_4_cleared.dot",
			"no_repeat_graph");

	if (output_contigs) {
		if (cfg::get().need_consensus) {
			string consensus_folder = cfg::get().output_dir
					+ "consensus_after_resolve/";
			OutputSingleFileContigs(resolved_gp.g, consensus_folder);
			string input_dir = cfg::get().input_dir;
			string reads_filename1 = input_dir + cfg::get().ds.first;
			string reads_filename2 = input_dir + cfg::get().ds.second;

			string real_reads = cfg::get().uncorrected_reads;
			if (real_reads != "none") {
				reads_filename1 = input_dir + cfg::get().ds.first;
				reads_filename2 = input_dir + cfg::get().ds.second;
			}

			typedef io::EasyReader EasyStream;
			EasyStream reads_1(reads_filename1);
			EasyStream reads_2(reads_filename2);
//			CarefulFilteringStream freads_1(reads_1);
//			CarefulFilteringStream freads_2(reads_2);
//			RCStream  frc_1(freads_1);
//			RCStream  frc_2(freads_2);
			vector<ReadStream*> reads = {/*&frc_1, &frc_2*/&reads_1, &reads_2 };

//			SelectReadsForConsensus<K,  graph_pack>(origin_gp, resolved_gp.g, labels_after, origin_gp.index, reads, consensus_folder);
			consensus_folder = cfg::get().output_dir
					+ "consensus_before_resolve/";
			OutputSingleFileContigs(origin_gp.g, consensus_folder);
			SelectReadsForConsensusBefore<K,  graph_pack>(origin_gp, origin_gp.g, labels_after, origin_gp.index, reads, consensus_folder);

		}

/*		one_many_contigs_enlarger<typename graph_pack::graph_t> N50enlarger(
				resolved_gp.g, cfg::get().ds.IS);
		N50enlarger.Loops_resolve();

		omnigraph::WriteSimple(
				resolved_gp.g,
				tot_labeler_after,
				cfg::get().output_dir + subfolder + graph_name
						+ "_5_unlooped.dot", "no_repeat_graph");

		OutputContigs(resolved_gp.g,
				cfg::get().output_dir + "unlooped" + postfix);
*/
//		N50enlarger.one_many_resolve_with_vertex_split();
//
//		omnigraph::WriteSimple(
//				resolved_gp.g,
//				tot_labeler_after,
//				cfg::get().output_dir + subfolder + graph_name
//						+ "_6_finished.dot", "no_repeat_graph");
//
//		OutputContigs(resolved_gp.g,
//				cfg::get().output_dir + "contigs_final.fasta");


	}
}

template<class graph_pack>
void component_statistics(graph_pack & conj_gp, int component_id,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index) {

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
	string component_name = cfg::get().output_dir + "graph_components/"
			+ graph_name;
	//component output
	string table_name = cfg::get().output_dir + "graph_components/tables/";
	make_dir(table_name);
	table_name += graph_name;
	set<typename graph_pack::graph_t::EdgeId> incoming_edges;
	set<typename graph_pack::graph_t::EdgeId> outgoing_edges;
	for (auto iter = conj_gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		typename graph_pack::graph_t::VertexId start = conj_gp.g.EdgeStart(
				*iter);
		typename graph_pack::graph_t::VertexId end = conj_gp.g.EdgeEnd(*iter);
		if (conj_gp.g.length(*iter) > cfg::get().ds.IS + 100) {

			if (conj_gp.g.IsDeadStart(
					start) /*&& conj_gp.g.CheckUniqueOutgoingEdge(start)*/) {
				incoming_edges.insert(*iter);
			} else if (conj_gp.g.IsDeadEnd(
					end)/* && conj_gp.g.CheckUniqueIncomingEdge(end)*/) {
				outgoing_edges.insert(*iter);
			} else {
				WARN(
						"strange long edge in component " << component_name << " , edge_id " << conj_gp.int_ids.ReturnIntId(*iter));
			}
		}
	}INFO("incoming- outgoint set formed");
	int flag = 1;
	for (auto inc_iter = incoming_edges.begin();
			inc_iter != incoming_edges.end(); ++inc_iter) {
		int count = 0;
		for (auto out_iter = outgoing_edges.begin();
				out_iter != outgoing_edges.end(); ++out_iter) {
			if (clustered_index.GetEdgePairInfo(*inc_iter, *out_iter).size()
					== 1)
				count++;
		}
		if (count != 1)
			flag = 0;
	}
	FILE* file;
	if (flag)
		file = fopen((table_name + ".tbl_good").c_str(), "w");
	else
		file = fopen((table_name + ".tbl").c_str(), "w");

	INFO("Saving in-out table , " << component_name <<" created");
	VERIFY(file != NULL);
	fprintf(file, "%7c", ' ');

	for (auto out_iter = outgoing_edges.begin();
			out_iter != outgoing_edges.end(); ++out_iter)
		fprintf(file, " %7d", conj_gp.int_ids.ReturnIntId(*out_iter));
	fprintf(file, "\n");

	for (auto inc_iter = incoming_edges.begin();
			inc_iter != incoming_edges.end(); ++inc_iter) {
		fprintf(file, " %7d", conj_gp.int_ids.ReturnIntId(*inc_iter));
		for (auto out_iter = outgoing_edges.begin();
				out_iter != outgoing_edges.end(); ++out_iter) {
			char c;
			if (clustered_index.GetEdgePairInfo(*inc_iter, *out_iter).size()
					== 0)
				c = '0';
			else
				c = 'X';
			fprintf(file, "%7c", c);
		}
		fprintf(file, "\n");
	}

	fprintf(file, "\n");
	for (auto inc_iter = incoming_edges.begin();
			inc_iter != incoming_edges.end(); ++inc_iter)
		fprintf(file, " %7d", conj_gp.int_ids.ReturnIntId(*inc_iter));
	fprintf(file, "\n");
	for (auto out_iter = outgoing_edges.begin();
			out_iter != outgoing_edges.end(); ++out_iter)
		fprintf(file, " %7d", conj_gp.int_ids.ReturnIntId(*out_iter));
	fprintf(file, "\n");

	fclose(file);

}

void resolve_conjugate_component(int component_id, const Sequence& genome) {
	conj_graph_pack conj_gp(genome);
	paired_info_index paired_index(conj_gp.g/*, 5.*/);
	paired_info_index clustered_index(conj_gp.g);

	INFO("Resolve component "<<component_id);

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
	string component_name = cfg::get().output_dir + "graph_components/"
			+ graph_name;

	ScanWithClusteredIndex(component_name, conj_gp, clustered_index);

	component_statistics(conj_gp, component_id, clustered_index);

	conj_graph_pack resolved_gp(genome);
	string sub_dir = "resolve_components/";

	string resolved_name = cfg::get().output_dir + "resolve_components"
			+ "/resolve_" + graph_name + "/";
	make_dir(resolved_name);
	process_resolve_repeats(conj_gp, clustered_index, resolved_gp, graph_name,
			sub_dir, false);
}

void resolve_nonconjugate_component(int component_id, const Sequence& genome) {
	nonconj_graph_pack nonconj_gp(genome);
	PairedInfoIndex<nonconj_graph_pack::graph_t> clustered_index(nonconj_gp.g);

	INFO("Resolve component "<<component_id);

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
	string component_name = cfg::get().output_dir + "graph_components/"
			+ graph_name;

	ScanWithClusteredIndex(component_name, nonconj_gp, clustered_index);

	component_statistics(nonconj_gp, component_id, clustered_index);

	nonconj_graph_pack resolved_gp(genome);
	string sub_dir = "resolve_components/";

	string resolved_name = cfg::get().output_dir + "resolve_components"
			+ "/resolve_" + graph_name + "/";
	make_dir(resolved_name);
	process_resolve_repeats(nonconj_gp, clustered_index, resolved_gp,
			graph_name, sub_dir, false);
}

Sequence load_genome() {
	string genome_filename = cfg::get().ds.reference_genome;
	std::string genome;
	if (genome_filename.length() > 0) {
		genome_filename = cfg::get().input_dir + genome_filename;
		checkFileExistenceFATAL(genome_filename);
		io::Reader<io::SingleRead> genome_stream(genome_filename);
		io::SingleRead full_genome;
		genome_stream >> full_genome;
		genome = full_genome.GetSequenceString().substr(0,
				cfg::get().ds.LEN); // cropped
	}
	return Sequence(genome);
}

void resolve_repeats() {
	// read data ('genome')
	Sequence genome = load_genome();

	conj_graph_pack conj_gp(genome);
	paired_info_index paired_index(conj_gp.g/*, 5.*/);
	paired_info_index clustered_index(conj_gp.g);

	exec_distance_estimation(conj_gp, paired_index, clustered_index);

	if(cfg::get().pos.late_threading){
		FillEdgesPos(conj_gp, conj_gp.genome, 10);
		FillEdgesPos(conj_gp, !conj_gp.genome, 11);
		FillEdgesPos(conj_gp, cfg::get().pos.contigs_for_threading, 10000);
	}
	
	//tSeparatedStats(conj_gp, conj_gp.genome, clustered_index);
	
	INFO("STAGE == Resolving Repeats");

	//todo refactor labeler creation
	total_labeler_graph_struct graph_struct(conj_gp.g, &conj_gp.int_ids, &conj_gp.edge_pos);
	total_labeler tot_lab(&graph_struct);
	EdgeQuality<Graph> quality_labeler(conj_gp.g, conj_gp.index, conj_gp.kmer_mapper, conj_gp.genome);
//	OutputWrongContigs<K>(conj_gp, 1000, "contamination.fasta");
	LabelerList<Graph> labeler(tot_lab, quality_labeler);
    detail_info_printer printer(conj_gp, labeler, cfg::get().output_dir, "graph.dot");
    printer(ipp_before_repeat_resolution);

	if (!cfg::get().paired_mode) {
		OutputContigs(conj_gp.g, cfg::get().output_dir + "contigs.fasta");
		return;
	}

	int number_of_components = 0;

	if (cfg::get().componential_resolve) {
		make_dir(cfg::get().output_dir + "graph_components" + "/");
		number_of_components = PrintGraphComponents(
				cfg::get().output_dir + "graph_components/graph_", conj_gp,
				cfg::get().ds.IS + 100, clustered_index);
		INFO("number of components "<<number_of_components);
	}

	if (cfg::get().rr.symmetric_resolve) {
		conj_graph_pack resolved_gp(genome);
		if (cfg::get().etalon_info_mode) {
			//temporary
			process_resolve_repeats(conj_gp, conj_gp.etalon_paired_index,
					resolved_gp, "graph");
		} else {
			process_resolve_repeats(conj_gp, clustered_index, resolved_gp,
					"graph");
		}
		if (cfg::get().componential_resolve) {
			make_dir(cfg::get().output_dir + "resolve_components" + "/");
			for (int i = 0; i < number_of_components; i++) {
				resolve_conjugate_component(i + 1, genome);
			}
		}
	} else {
		nonconj_graph_pack origin_gp(conj_gp.genome);
		PairedInfoIndex<nonconj_graph_pack::graph_t> orig_clustered_idx(origin_gp.g);
		Convert(conj_gp, clustered_index, origin_gp, orig_clustered_idx);
		nonconj_graph_pack resolved_gp(conj_gp.genome);
		process_resolve_repeats(origin_gp, orig_clustered_idx,
				resolved_gp, "graph");
		if (cfg::get().componential_resolve) {
			make_dir(cfg::get().output_dir + "resolve_components" + "/");
			for (int i = 0; i < number_of_components; i++) {
				resolve_nonconjugate_component(i + 1, genome);
			}
		}
	}

}

void exec_repeat_resolving() {
	if (cfg::get().entry_point <= ws_repeats_resolving) {
		resolve_repeats();
		//todo why nothing to save???
		// nothing to save yet
	} else {
		INFO("Loading Repeat Resolving");
		INFO("Nothing to load");
		// nothing to load
	}
}

} // debruijn_graph

