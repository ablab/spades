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

class PacBioAligner {
public:

private:

	conj_graph_pack &gp_;
	size_t k_test_;DECL_LOGGER("PacIndex")

public:
	typedef Graph::EdgeId EdgeId;

	PacBioAligner(conj_graph_pack& conj_gp, size_t k_test) :
			gp_(conj_gp), k_test_(k_test) {
	}

	void pacbio_test(PathStorage<Graph> &long_reads, GapStorage<Graph> &gaps) {
		INFO("starting pacbio tests");
		ReadStream* pacbio_read_stream = new io::EasyReader(cfg::get().pacbio_reads, false);
		size_t n = 0;
		//    map<int, int> profile;
		map<int, int> different_edges_profile;
		int genomic_subreads = 0;
		int nongenomic_subreads = 0;
		int nongenomic_edges = 0;
		int total_length = 0;
		int tlen = 0;
//	    PathStorage<Graph> long_reads(gp_.g);
//	    long_reads.LoadFromFile("long_reads.mpr");
//	    INFO("dumping back");
//	    long_reads.DumpToFile("long_reads2.mpr", gp_.edge_pos);

		int rc_pairs = 0;
		size_t read_buffer_size = 50000;
		std::vector<ReadStream::read_type> reads(read_buffer_size);
		ReadStream::read_type read;
		size_t buffer_no = 0;
		PacBioMappingIndex<ConjugateDeBruijnGraph> pac_index(gp_.g, k_test_);
		ofstream filestr("pacbio_mapped.mpr");
		filestr.close();
		while (!pacbio_read_stream->eof()) {
			size_t buf_size = 0;
			for (; buf_size < read_buffer_size && !pacbio_read_stream->eof(); ++buf_size) {
				*pacbio_read_stream >> reads[buf_size];
			}
			INFO("Prepared batch " << buffer_no << " of " << buf_size << " reads.");
			DEBUG("master thread number " << omp_get_thread_num());
			ProcessReadsBatch(reads, pac_index, long_reads, gaps, buf_size, genomic_subreads, nongenomic_subreads, nongenomic_edges, total_length, tlen, n,
					different_edges_profile, rc_pairs);
			INFO("Processed batch " << buffer_no);
			++buffer_no;
		}

		long_reads.DumpToFile("long_reads.mpr", gp_.edge_pos);
		gaps.DumpToFile("gaps.mpr", gp_.edge_pos);
		gaps.PadGapStrings();
		gaps.DumpToFile("gaps_padded.mpr", gp_.edge_pos);
		PacbioGapCloser<Graph> gap_closer(gp_.g);
		gap_closer.ConstructConsensus(cfg::get().max_threads, gaps);
		gap_closer.DumpToFile("gaps_closed.fasta", gp_.edge_pos);
		INFO("Total reads: " << n);
		INFO("Mean read length: " << total_length * 0.1/ n)
		INFO("Mean subread length: " << tlen * 0.1/ (genomic_subreads + nongenomic_subreads + nongenomic_edges))
		INFO("reads with rc edges:  " << rc_pairs);
		INFO("Genomic/nongenomic subreads/nongenomic edges: "<<genomic_subreads <<" / " << nongenomic_subreads <<" / "<< nongenomic_edges);
		INFO("different edges profile:")
		for (auto iter = different_edges_profile.begin(); iter != different_edges_profile.end(); ++iter)
			if (iter->first < 100) {
				INFO(iter->first <<" :  "<< iter->second);
			}
		INFO("PacBio test finished");
//		return ;
	}

	void ProcessReadsBatch(std::vector<ReadStream::read_type>& reads, PacBioMappingIndex<ConjugateDeBruijnGraph>& pac_index, PathStorage<Graph>& long_reads, GapStorage<Graph>& gaps,
			size_t buf_size, int& genomic_subreads, int& nongenomic_subreads, int& nongenomic_edges, int& total_length, int& tlen, size_t& n,
			map<int, int>& different_edges_profile, int& rc_pairs) {
		omp_lock_t tmp_file_output;
		omp_init_lock(&tmp_file_output);
		ofstream filestr("pacbio_mapped.mpr", ofstream::app);
		vector <PathStorage<Graph > > long_reads_by_thread(cfg::get().max_threads, PathStorage<Graph>(gp_.g));
		vector <GapStorage<Graph > > gaps_by_thread(cfg::get().max_threads, GapStorage<Graph>(gp_.g));

# pragma omp parallel for shared(reads, long_reads_by_thread, pac_index, n,  different_edges_profile) num_threads(cfg::get().max_threads)
		for (size_t i = 0; i < buf_size; ++i) {
			if (i % 1000 == 0) {
				DEBUG("thread number " << omp_get_thread_num());
			}
			size_t thread_num = omp_get_thread_num();
			Sequence seq(reads[i].sequence());
			total_length += seq.size();
			auto location_map = pac_index.GetClusters(seq);
			different_edges_profile[location_map.size()]++;
			n++;
//		    if (location_map.size() <= 1){
//		    	TRACE("No significant clusters");
//		    	continue;
//		    }
			for (auto iter = location_map.begin(); iter != location_map.end(); ++iter) {
				bool flag = false;
				for (auto j_iter = location_map.begin(); j_iter != location_map.end(); ++j_iter) {
					if (iter != j_iter && gp_.g.conjugate(iter->edgeId) == j_iter->edgeId) {
						flag = true;
						break;
					}
				}
				if (flag == true) {
					rc_pairs++;
					break;
				}
			}
			DEBUG(n << "  " << location_map.size()<< ": \n");
			auto current_read_mapping = pac_index.GetReadAlignment(seq);
			auto aligned_edges = current_read_mapping.main_storage;
			auto gaps = current_read_mapping.gaps;
			for (auto iter = gaps.begin(); iter != gaps.end(); ++iter) {
				gaps_by_thread[thread_num].AddGap(*iter, true);
			}
			for (auto iter = aligned_edges.begin(); iter != aligned_edges.end(); ++iter) {
				long_reads_by_thread[thread_num].AddPath(*iter, 1, true);
				if (gp_.edge_pos.IsConsistentWithGenome(*iter)) {
					genomic_subreads++;
				} else {
					if (iter->size() > 1)
						nongenomic_subreads++;
					else
						nongenomic_edges++;
				}
			}
//this block is something to be overcome
			omp_set_lock(&tmp_file_output);
			filestr << n << "  " << location_map.size() << ": \n";
			for (auto iter = location_map.begin(); iter != location_map.end(); ++iter) {
				filestr << gp_.g.int_id(iter->edgeId) << "(" << gp_.g.length(iter->edgeId) << ")  " << iter->sorted_positions.size() << "\n";
				for (auto set_iter = iter->sorted_positions.begin(); set_iter != iter->sorted_positions.end(); ++set_iter)
					filestr << set_iter->edge_position << "-" << set_iter->read_position << "   ";
				filestr << " \n";
			}
			filestr << "found " << aligned_edges.size() << " aligned subreads.\n";
			for (auto iter = aligned_edges.begin(); iter != aligned_edges.end(); ++iter) {
				string tmp = " ";
				if (gp_.edge_pos.IsConsistentWithGenome(*iter)) {

				} else {
					tmp = " NOT ";
				}
				filestr << "Alignment of " << iter->size() << " edges is" << tmp << "consistent with genome\n";
				//except this point
				for (auto j_iter = iter->begin(); j_iter != iter->end(); ++j_iter) {
					filestr << gp_.g.int_id(*j_iter) << "(" << gp_.g.length(*j_iter) << ") ";
					tlen += gp_.g.length(*j_iter);
				}
				filestr << " \n";
			}
			filestr << " \n";
			filestr << " \n";
			omp_unset_lock(&tmp_file_output);

			VERBOSE_POWER(n, " reads processed");

		}
		for (size_t i = 0; i < cfg::get().max_threads; i++) {
			long_reads.AddStorage(long_reads_by_thread[i]);
			gaps.AddStorage(gaps_by_thread[i]);
		}
//			filestr << "Gap: " << gp_.g.int_id(iter->first.first) << " " << gp_.g.int_id(iter->first.second) << " weight  " << iter->second << "\n";

		filestr.close();
	}

};

