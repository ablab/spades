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
#include "io/careful_filtering_reader_wrapper.hpp"
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

namespace debruijn_graph
{

void FillContigNumbers(   map<ConjugateDeBruijnGraph::EdgeId, int>& contigNumbers,  ConjugateDeBruijnGraph& cur_graph){
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

void FillContigNumbers(   map<NonconjugateDeBruijnGraph::EdgeId, int>& contigNumbers,  NonconjugateDeBruijnGraph& cur_graph){
	int cur_num = 0;
	for (auto iter = cur_graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		contigNumbers[*iter] = cur_num;
	    cur_num++;
	}
}

template<size_t k, class Graph>
void SelectReadsForConsensus(Graph& etalon_graph, Graph& cur_graph,
        EdgeLabelHandler<Graph>& LabelsAfter,
        const EdgeIndex<K + 1, Graph>& index ,vector<ReadStream *>& reads
        , string& consensus_output_dir)
{
    INFO("ReadMapping started");
    map<typename Graph::EdgeId, int> contigNumbers;
    int cur_num = 0;
    FillContigNumbers(contigNumbers, cur_graph);
//    for (auto iter = cur_graph.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
//        contigNumbers[*iter] = cur_num;
//        cur_num++;
//    }
//
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
        SingleReadMapper<k, Graph> rm(etalon_graph, index);
        INFO("mapping reads from pair"<< i);
        while (!reads[i - 1]->eof()) {
            io::SingleRead cur_read;

            (* reads[i - 1]) >> cur_read;
            vector<typename Graph::EdgeId> res = rm.GetContainingEdges(
                    cur_read);
            read_num++;
            TRACE(
                    read_num<< " mapped to"<< res.size() <<" contigs :, read"<< cur_read.sequence());
//          map_quantity += res.size();
            for (size_t ii = 0; ii < res.size(); ii++) {
                TRACE("conting number "<< contigNumbers[res[ii]]);
                set<typename Graph::EdgeId> images =
                        LabelsAfter.edge_inclusions[res[ii]];
                for (auto iter = images.begin(); iter != images.end(); ++iter)
                    (*mapped_reads[contigNumbers[*iter]])
                            << cur_read.sequence();
            }
        }
    }
}


template <class graph_pack>
void process_resolve_repeats(
		graph_pack& origin_gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
		graph_pack& resolved_gp, const string& graph_name, const string& subfolder = "", bool output_contigs = true)
{
    EdgeLabelHandler<typename graph_pack::graph_t> labels_after(resolved_gp.g, origin_gp.g);
    DEBUG("New index size: "<< clustered_index.size());
    // todo: make printGraph const to its arguments

    // todo: possibly we don't need it
//    if (cfg::get().rectangle_mode)
//        RectangleResolve(clustered_index, origin_gp.g, cfg::get().output_root + "tmp/", cfg::get().output_dir);

    typedef TotalLabelerGraphStruct<typename graph_pack::graph_t> total_labeler_gs;
    typedef TotalLabeler           <typename graph_pack::graph_t> total_labeler;


    total_labeler_gs graph_struct_before(origin_gp  .g, &origin_gp  .int_ids, &origin_gp  .edge_pos, NULL);
    total_labeler tot_labeler_before(&graph_struct_before);

    omnigraph::WriteSimple(origin_gp.g, tot_labeler_before, cfg::get().output_dir + subfolder + graph_name + "_2_simplified.dot", "no_repeat_graph");

    ResolveRepeats(origin_gp  .g, origin_gp  .int_ids, clustered_index, origin_gp  .edge_pos,
                   resolved_gp.g, resolved_gp.int_ids,                  resolved_gp.edge_pos,
                   cfg::get().output_dir + subfolder +"resolve_" + graph_name +  "/", labels_after);
    if (output_contigs) {
       	OutputContigs(resolved_gp.g, cfg::get().output_dir + "contigs_after_rr_before_simplify.fasta");
    	OutputContigs(origin_gp.g, cfg::get().output_dir + "contigs_before_resolve.fasta");
    }
    INFO("Total labeler start");

    total_labeler_gs graph_struct_after (resolved_gp.g, &resolved_gp.int_ids, &resolved_gp.edge_pos, &labels_after);

    total_labeler tot_labeler_after(&graph_struct_after, &graph_struct_before);

    omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after, cfg::get().output_dir + subfolder + graph_name + "_3_resolved.dot", "no_repeat_graph");

    INFO("Total labeler finished");

    INFO("---Clearing resolved graph---");


    for (int i = 0; i < 3; ++i)
    {
        ClipTipsForResolver(resolved_gp.g);
        BulgeRemoveWrap      (resolved_gp.g);
        RemoveLowCoverageEdges(resolved_gp.g, i, 3);
//        RemoveRelativelyLowCoverageEdges(resolved_gp.g);
    }

    INFO("---Cleared---");
    INFO("---Output Contigs---");

    if (output_contigs)
    	OutputContigs(resolved_gp.g, cfg::get().output_dir + "contigs_before_enlarge.fasta");

    omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,


    cfg::get().output_dir + subfolder + graph_name + "_4_cleared.dot", "no_repeat_graph");

    if (output_contigs)
    {
		if (cfg::get().need_consensus)
		{
			string consensus_folder = cfg::get().output_dir + "consensus_after_resolve/";
			OutputSingleFileContigs(resolved_gp.g, consensus_folder);
			string input_dir = cfg::get().input_dir;
			string reads_filename1 = input_dir + cfg::get().ds.first;
			string reads_filename2 = input_dir + cfg::get().ds.second;

			string real_reads = cfg::get().uncorrected_reads;
			if (real_reads != "none") {
				reads_filename1 = input_dir + cfg::get().ds.first;
				reads_filename2 = input_dir + cfg::get().ds.second;
			}

			typedef io::EasyReader<io::SingleRead> EasyStream;
			EasyStream reads_1(reads_filename1);
			EasyStream reads_2(reads_filename2);
//			CarefulFilteringStream freads_1(reads_1);
//			CarefulFilteringStream freads_2(reads_2);
//			RCStream  frc_1(freads_1);
//			RCStream  frc_2(freads_2);
			vector<ReadStream*> reads = {/*&frc_1, &frc_2*/&reads_1, &reads_2};

			SelectReadsForConsensus<K, typename graph_pack::graph_t>(origin_gp.g, resolved_gp.g, labels_after, origin_gp.index, reads, consensus_folder);
		}

		one_many_contigs_enlarger<typename graph_pack::graph_t> N50enlarger(resolved_gp.g, cfg::get().ds.IS);
		N50enlarger.Loops_resolve();

		omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,
							 cfg::get().output_dir  + subfolder + graph_name + "_5_unlooped.dot", "no_repeat_graph");

		OutputContigs(resolved_gp.g, cfg::get().output_dir + "contigs_unlooped.fasta");

		N50enlarger.one_many_resolve_with_vertex_split();

		omnigraph::WriteSimple(resolved_gp.g, tot_labeler_after,
							 cfg::get().output_dir + subfolder + graph_name + "_6_finished.dot", "no_repeat_graph");

		OutputContigs(resolved_gp.g, cfg::get().output_dir + "contigs_final.fasta");
		OutputContigs(origin_gp.g, cfg::get().output_dir + "contigs_before_resolve.fasta");

    }
}


template <class graph_pack>
void component_statistics(graph_pack & conj_gp, int component_id, PairedInfoIndex<typename graph_pack::graph_t>& clustered_index){

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
	string component_name = cfg::get().output_dir + "graph_components/" + graph_name;
	//component output
	string table_name = cfg::get().output_dir + "graph_components/tables/";
	mkdir(table_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	table_name += graph_name;
	set<typename graph_pack::graph_t::EdgeId> incoming_edges;
	set<typename graph_pack::graph_t::EdgeId> outgoing_edges;
	for(auto iter = conj_gp.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		typename graph_pack::graph_t::VertexId start = conj_gp.g.EdgeStart(*iter);
		typename graph_pack::graph_t::VertexId end = conj_gp.g.EdgeEnd(*iter);
		if (conj_gp.g.length(*iter) > cfg::get().ds.IS + 100) {

			if (conj_gp.g.IsDeadStart(start) /*&& conj_gp.g.CheckUniqueOutgoingEdge(start)*/){
				incoming_edges.insert(*iter);
			} else if (conj_gp.g.IsDeadEnd(end)/* && conj_gp.g.CheckUniqueIncomingEdge(end)*/){
				outgoing_edges.insert(*iter);
			} else {
				WARN("strange long edge in component " << component_name << " , edge_id " << conj_gp.int_ids.ReturnIntId(*iter));
			}
		}
	}
	INFO("incoming- outgoint set formed");
	int flag = 1;
	for (auto inc_iter = incoming_edges.begin(); inc_iter != incoming_edges.end(); ++inc_iter){
		int count = 0;
 		for (auto out_iter = outgoing_edges.begin(); out_iter != outgoing_edges.end(); ++out_iter){
			if (clustered_index.GetEdgePairInfo(*inc_iter, *out_iter).size() == 1)
				count ++;
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
	fprintf(file,"%7c", ' ');

	for (auto out_iter = outgoing_edges.begin(); out_iter != outgoing_edges.end(); ++out_iter)
		fprintf(file," %7d", conj_gp.int_ids.ReturnIntId(*out_iter));
	fprintf(file, "\n");

	for (auto inc_iter = incoming_edges.begin(); inc_iter != incoming_edges.end(); ++inc_iter){
		fprintf(file," %7d", conj_gp.int_ids.ReturnIntId(*inc_iter));
		for (auto out_iter = outgoing_edges.begin(); out_iter != outgoing_edges.end(); ++out_iter){
			char c;
			if (clustered_index.GetEdgePairInfo(*inc_iter, *out_iter).size() == 0)
				c = '0';
			else
				c = 'X';
			fprintf(file,"%7c", c);
		}
		fprintf(file, "\n");
	}

	fprintf(file, "\n");
	for (auto inc_iter = incoming_edges.begin(); inc_iter != incoming_edges.end(); ++inc_iter)
		fprintf(file," %7d", conj_gp.int_ids.ReturnIntId(*inc_iter));
	fprintf(file, "\n");
	for (auto out_iter = outgoing_edges.begin(); out_iter != outgoing_edges.end(); ++out_iter)
		fprintf(file," %7d", conj_gp.int_ids.ReturnIntId(*out_iter));
	fprintf(file, "\n");

	fclose(file);

}



void resolve_conjugate_component(int component_id, const Sequence& genome){
    conj_graph_pack   conj_gp (genome);
    paired_info_index paired_index   (conj_gp.g/*, 5.*/);
    paired_info_index clustered_index(conj_gp.g);

    INFO("Resolve component "<<component_id);

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
	string component_name = cfg::get().output_dir + "graph_components/" + graph_name;

	scanConjugateGraph(&conj_gp.g, &conj_gp.int_ids, component_name, &clustered_index,
			&conj_gp.edge_pos);

	component_statistics(conj_gp, component_id, clustered_index);


	conj_graph_pack   resolved_gp (genome);
    string sub_dir = "resolve_components/";


    string resolved_name = cfg::get().output_dir + "resolve_components" + "/resolve_" + graph_name + "/";
	mkdir(resolved_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	process_resolve_repeats(conj_gp, clustered_index, resolved_gp, graph_name, sub_dir, false ) ;
}

void resolve_nonconjugate_component(int component_id, const Sequence& genome){
    nonconj_graph_pack   nonconj_gp;
//    paired_info_index paired_index   (nonconj_gp.g/*, 5.*/);

    INFO("Resolve component "<<component_id);

	string graph_name = ConstructComponentName("graph_", component_id).c_str();
	string component_name = cfg::get().output_dir + "graph_components/" + graph_name;

	scanNCGraph(nonconj_gp.g, nonconj_gp.int_ids, component_name, &nonconj_gp.clustered_index,
			nonconj_gp.edge_pos);

	component_statistics(nonconj_gp, component_id, nonconj_gp.clustered_index);


	nonconj_graph_pack   resolved_gp;
	string sub_dir = "resolve_components/";


	string resolved_name = cfg::get().output_dir + "resolve_components" + "/resolve_" + graph_name + "/";
	mkdir(resolved_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	process_resolve_repeats(nonconj_gp, nonconj_gp.clustered_index, resolved_gp, graph_name, sub_dir, false ) ;
}




void resolve_repeats(PairedReadStream& stream, const Sequence& genome)
{
    conj_graph_pack   conj_gp (genome);
    paired_info_index paired_index   (conj_gp.g/*, 5.*/);
    paired_info_index clustered_index(conj_gp.g);

    exec_distance_estimation(stream, conj_gp, paired_index, clustered_index);

    INFO("STAGE == Resolving Repeats");

    if (!cfg::get().paired_mode)
    {
        OutputContigs(conj_gp.g, cfg::get().output_dir + "contigs.fasta");
        return;
    }

    int number_of_components = 0;


    if (cfg::get().componential_resolve){
		mkdir((cfg::get().output_dir + "graph_components" + "/").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
		number_of_components = PrintGraphComponents(
				cfg::get().output_dir + "graph_components/graph_", conj_gp.g,
				cfg::get().ds.IS+100, conj_gp.int_ids, clustered_index, conj_gp.edge_pos,
				cfg::get().rr.symmetric_resolve);
		INFO("number of components "<<number_of_components);
    }


    if (cfg::get().rr.symmetric_resolve) {
    	conj_graph_pack   resolved_gp (genome);
        if (cfg::get().etalon_info_mode){
        	//temporary
        	process_resolve_repeats(conj_gp, conj_gp.etalon_paired_index, resolved_gp, "graph") ;
        } else {
        	process_resolve_repeats(conj_gp, clustered_index, resolved_gp, "graph") ;
        }
        if (cfg::get().componential_resolve){
        	mkdir((cfg::get().output_dir + "resolve_components" + "/").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
        	for (int i = 0; i<number_of_components; i++){
        		resolve_conjugate_component(i+1, genome);
        	}
        }
    } else {
    	nonconj_graph_pack origin_gp(conj_gp, clustered_index);
    	nonconj_graph_pack resolved_gp;
       	process_resolve_repeats(origin_gp, origin_gp.clustered_index, resolved_gp, "graph") ;
        if (cfg::get().componential_resolve){
        	mkdir((cfg::get().output_dir + "resolve_components" + "/").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
        	for (int i = 0; i<number_of_components; i++){
        		resolve_nonconjugate_component(i+1, genome);
        	}
        }
    }

}


void exec_repeat_resolving(PairedReadStream& stream, const Sequence& genome)
{
	if (cfg::get().entry_point <= ws_repeats_resolving)
    {
		resolve_repeats(stream, genome);
		//todo why nothing to save???
        // nothing to save yet
    }
    else
    {
    	INFO("Loading Repeat Resolving");
    	INFO("Nothing to load");
    	// nothing to load
    }
}

} // debruijn_graph






