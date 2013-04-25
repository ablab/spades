/*
 * pacbio_aligner.hpp
 *
 *  Created on: Apr 25, 2013
 *      Author: lab42
 */

#pragma once

#include "debruijn_kmer_index.hpp"
#include "graph_pack.hpp"
#include <algorithm>
#include "pac_index.hpp"
#include "long_read_storage.hpp"

class PacBioAligner{
public:


protected :
	conj_graph_pack &gp_;
	size_t k_test_;
	PacBioMappingIndex<ConjugateDeBruijnGraph> pac_index;

private:
	DECL_LOGGER("PacIndex")

public :
	typedef typename Graph::EdgeId EdgeId;

	PacBioAligner(conj_graph_pack& conj_gp, size_t k_test):gp_(conj_gp), k_test_(k_test),pac_index(gp_.g, k_test){}


	void pacbio_test(){

	    INFO("starting pacbio tests");
		ofstream filestr("pacbio_mapped.mpr");

		ReadStream* pacbio_read_stream = new io::EasyReader(cfg::get().pacbio_reads, true);
	    size_t n = 0;
	//    map<int, int> profile;
	    map<int, int> different_edges_profile;
	    int genomic_subreads = 0;
	    int nongenomic_subreads = 0;
	    int nongenomic_edges = 0;
	    int total_length = 0;
	    int tlen = 0;
	    LongReadStorage<Graph> long_reads(gp_.g);
	    int rc_pairs = 0;
		while (!pacbio_read_stream->eof()) {
			ReadStream::read_type read;
			*pacbio_read_stream>>read;
			Sequence seq(read.sequence());
			total_length += seq.size();
		    auto location_map = pac_index.GetClusters(seq);
		    //	    size_t res_count = pac_index.Count(seq);
		    //	    if (profile.find(res_count) == profile.end())
		    //	    	profile.insert(make_pair(res_count, 0));
		    //	    profile[res_count] ++;
		    //	    if (res_count != 0){
		    //	    	DEBUG(read.sequence());
		    //	    	DEBUG(res_count);
		    //	    }

		    different_edges_profile[location_map.size()]++;
		    n++;
		    if (location_map.size() <= 1){
		    	TRACE("No significant clusters");
		    	continue;
		    }
		    for (auto iter = location_map.begin(); iter != location_map.end(); ++iter) {
		    	bool flag = false;
			    for (auto j_iter = location_map.begin(); j_iter != location_map.end(); ++j_iter) {
			    	if (iter != j_iter && gp_.g.conjugate(iter->first) == j_iter->first) {
			    		flag = true;
			    		break;
			    	}
			    }
		    	if (flag == true) {
		    		rc_pairs ++;
		    		break;
		    	}
		    }
		    //continue;
		    filestr << n << "  " << location_map.size()<< ": \n";
		    INFO(n << "  " << location_map.size()<< ": \n");
		    for (auto iter = location_map.begin(); iter != location_map.end(); ++iter) {
		    	filestr << gp_.g.int_id(iter->first) <<"("<<gp_.g.length(iter->first)<<")  " << iter->second.size() <<"\n";
		    	for (auto set_iter = iter->second.begin(); set_iter != iter->second.end(); ++ set_iter)
					filestr << set_iter->edge_position << "-" << set_iter->read_position << "   ";
				filestr << " \n";
		    }

		    auto aligned_edges = pac_index.GetReadAlignment(seq);
		    filestr <<"found "<< aligned_edges.size()  <<" aligned subreads.\n";
		    for(auto iter = aligned_edges.begin(); iter != aligned_edges.end(); ++iter) {
		    	string tmp = " ";
		    	if (gp_.edge_pos.IsConsistentWithGenome(*iter)) {
		    		genomic_subreads ++;
		    	}else {
		    		tmp = " NOT ";
		    		if (iter->size() > 1)
		    			nongenomic_subreads ++;
		    		else
		    			nongenomic_edges ++;
		    	}
		    	filestr <<"Alignment of "<< iter->size()  <<" edges is" << tmp <<"consistent with genome\n";
		    	long_reads.AddPath(*iter);
		    	for (auto j_iter = iter->begin(); j_iter != iter->end(); ++j_iter){
		    		filestr << gp_.g.int_id(*j_iter) <<"("<<gp_.g.length(*j_iter)<<") ";
		    		tlen += gp_.g.length(*j_iter);
		    	}
		    	filestr << " \n";
		    }

		    filestr << " \n";
		    filestr << " \n";
		    VERBOSE_POWER(n, " reads processed");

		}
		long_reads.DumpToFile("long_reads.mpr", gp_.edge_pos);
		INFO("Total reads: " << n);
		INFO("Mean read length: " << total_length * 0.1/ n)
		INFO("Mean subread length: " << tlen * 0.1/ (genomic_subreads + nongenomic_subreads + nongenomic_edges))
		INFO("reads with rc edges:  " << rc_pairs);
		INFO("Genomic/nongenomic subreads/nongenomic edges: "<<genomic_subreads <<" / " << nongenomic_subreads <<" / "<< nongenomic_edges);
	//	INFO("profile:")
	//	for (auto iter = profile.begin(); iter != profile.end(); ++iter)
	//		if (iter->first < 100) {
	//			INFO(iter->first <<" :  "<< iter->second);
	//		}

		INFO("different edges profile:")
		for (auto iter = different_edges_profile.begin(); iter != different_edges_profile.end(); ++iter)
			if (iter->first < 100) {
				INFO(iter->first <<" :  "<< iter->second);
			}


		INFO("PacBio test finished");

	}

};



