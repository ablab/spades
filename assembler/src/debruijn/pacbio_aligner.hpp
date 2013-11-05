/*
 * pacbio_aligner.hpp
 *
 *  Created on: Apr 25, 2013
 *      Author: lab42
 */

#pragma once

#include <cstdlib>
#include "indices/debruijn_kmer_index.hpp"
#include "graph_pack.hpp"
#include <algorithm>
#include "pac_index.hpp"
#include "pacbio_gap_closer.hpp"
#include "long_read_storage.hpp"
#include "path_extend/pe_io.hpp"

class PacBioAligner {
public:

private:

	conj_graph_pack &gp_;
	PairedIndicesT& paired_indices;
	PairedIndicesT& clustered_indices;
	PairedIndicesT& scaffold_indices_;
	LongReadContainerT& single_long_reads_;
	size_t k_test_;DECL_LOGGER("PacIndex")

public:
	typedef Graph::EdgeId EdgeId;

	PacBioAligner(conj_graph_pack& conj_gp, PairedIndicesT& paired_indices,
                  PairedIndicesT& clustered_indices,
                  PairedIndicesT& scaffold_indices,
                  LongReadContainerT& single_long_reads,
                  size_t k_test)
            : gp_(conj_gp),
              paired_indices(paired_indices),
              clustered_indices(clustered_indices),
              scaffold_indices_(scaffold_indices),
              single_long_reads_(single_long_reads),
              k_test_(k_test) {
    }

	void pacbio_test(PathStorage<Graph> &long_reads, GapStorage<Graph> &gaps) {
	    if (cfg::get().entry_point <= ws_pacbio_aligning) {
            INFO("starting pacbio tests");
            ReadStream* pacbio_read_stream = new io::EasyReader(cfg::get().pb.pacbio_reads, false);
            size_t n = 0;
            //    map<int, int> profile;
            map<int, int> different_edges_profile;
            int genomic_subreads = 0;
            int nongenomic_subreads = 0;
            int nongenomic_edges = 0;
            int total_length = 0;
            int tlen = 0;
        //		gaps.LoadFromFile("gaps_padded.mpr");
        //		gaps.PostProcess();
        //		PacbioGapCloser<Graph> gap1_closer(gp_.g);
        //		gap1_closer.ConstructConsensus(cfg::get().max_threads, gaps);
        //		gap1_closer.DumpToFile("gaps_closed2.fasta", gp_.edge_pos);
        //		INFO("PacBio test finished");
        //		exit(0);
            int rc_pairs = 0;
            size_t read_buffer_size = 50000;
            std::vector<ReadStream::read_type> reads(read_buffer_size);
            ReadStream::read_type read;
            size_t buffer_no = 0;
            PacBioMappingIndex<ConjugateDeBruijnGraph> pac_index(gp_.g, k_test_, cfg::get().K);

            path_extend::ContigWriter cw(gp_.g);
            cw.writeEdges("before_rr_with_ids.fasta");
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
            map<EdgeId, EdgeId> replacement;
            long_reads.DumpToFile(cfg::get().output_saves +  "long_reads_before_rep.mpr", gp_.edge_pos, replacement);
            gaps.DumpToFile(cfg::get().output_saves + "gaps.mpr");
            gaps.PadGapStrings();
            gaps.DumpToFile(cfg::get().output_saves +  "gaps_padded.mpr");
            PacbioGapCloser<Graph> gap_closer(gp_.g);
            gap_closer.ConstructConsensus(cfg::get().max_threads, gaps);
            gap_closer.CloseGapsInGraph(replacement);
            long_reads.ReplaceEdges(replacement);
        //		gp_.edge_pos.clear();
        //		FillPos(gp_, gp_.genome, "10");
        //		FillPos(gp_, !gp_.genome, "11");

            //long_reads.DumpToFile(cfg::get().output_saves + "pacbio_aligned.mpr", gp_.edge_pos, replacement);
            save_pacbio_aligned(long_reads, replacement);
            gap_closer.DumpToFile(cfg::get().output_saves+ "gaps_pb_closed.fasta", gp_.edge_pos);
            INFO("Total reads: " << n);
            INFO("Mean read length: " << total_length * 0.1 / (double) n)
            INFO("Mean subread length: " << tlen * 0.1 / (genomic_subreads + nongenomic_subreads + nongenomic_edges))
            INFO("reads with rc edges:  " << rc_pairs);
            INFO("Genomic/nongenomic subreads/nongenomic edges: "<<genomic_subreads <<" / " << nongenomic_subreads <<" / "<< nongenomic_edges);
            INFO("different edges profile:")
            for (auto iter = different_edges_profile.begin(); iter != different_edges_profile.end(); ++iter)
                if (iter->first < 100) {
                    INFO(iter->first <<" :  "<< iter->second);
                }

            INFO("PacBio test finished");
            return ;
	    } else {
//	        INFO("Loading pacbio_aligned");
//            path::files_t used_files;
//            load_pacbio_aligned(long_reads, &used_files);
//            link_files_by_prefix(used_files, cfg::get().output_saves);
	    }
	}

	void save_pacbio_aligned(const PathStorage<Graph>& long_reads, map<EdgeId, EdgeId>& replacement)
	{
	  if (cfg::get().make_saves) {
	    string p = path::append_path(cfg::get().output_saves, "pacbio_aligning");
	    INFO("Saving current state to " << p);
        long_reads.DumpToFile(p + ".mpr", gp_.edge_pos, replacement);
	    PrintAll(p, gp_, paired_indices, clustered_indices, scaffold_indices_, single_long_reads_);
	    write_lib_data(p);
	  }
	}

	void load_pacbio_aligned(PathStorage<Graph>& long_reads, path::files_t* used_files)
	{
	    string p = path::append_path(cfg::get().load_from, "pacbio_aligning");
	    used_files->push_back(p);
	    ScanAll(p, gp_, paired_indices, clustered_indices, scaffold_indices_, single_long_reads_);
	    long_reads.LoadFromFile(p + ".mpr");
	    load_lib_data(p);
	}


	void ProcessReadsBatch(std::vector<ReadStream::read_type>& reads, PacBioMappingIndex<ConjugateDeBruijnGraph>& pac_index, PathStorage<Graph>& long_reads, GapStorage<Graph>& gaps,
			size_t buf_size, int& genomic_subreads, int& nongenomic_subreads, int& nongenomic_edges, int& total_length, int& tlen, size_t& n,
			map<int, int>& different_edges_profile, int& rc_pairs) {
		omp_lock_t tmp_file_output;
		omp_init_lock(&tmp_file_output);
		ofstream filestr("pacbio_mapped.mpr", ofstream::app);
		vector <PathStorage<Graph > > long_reads_by_thread(cfg::get().max_threads, PathStorage<Graph>(gp_.g));
		vector <GapStorage<Graph > > gaps_by_thread(cfg::get().max_threads, GapStorage<Graph>(gp_.g));
		GenomeConsistenceChecker<typename GraphPack::graph_t> checker(gp, 10, 0.2);

		# pragma omp parallel for shared(reads, long_reads_by_thread, pac_index, n,  different_edges_profile) num_threads(cfg::get().max_threads)
		for (size_t i = 0; i < buf_size; ++i) {
			if (i % 1000 == 0) {
				DEBUG("thread number " << omp_get_thread_num());
			}
			size_t thread_num = omp_get_thread_num();
			Sequence seq(reads[i].sequence());
			total_length += (int) seq.size();
			n++;
			auto location_map = pac_index.GetClusters(seq);
			different_edges_profile[(int) location_map.size()]++;
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
				if (checker.IsConsistentWithGenome(*iter)) {
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
				if (checker.IsConsistentWithGenome(*iter)) {

				} else {
					tmp = " NOT ";
				}
				filestr << "Alignment of " << iter->size() << " edges is" << tmp << "consistent with genome\n";
				//except this point
				for (auto j_iter = iter->begin(); j_iter != iter->end(); ++j_iter) {
					filestr << gp_.g.int_id(*j_iter) << "(" << gp_.g.length(*j_iter) << ") ";
					tlen += (int) gp_.g.length(*j_iter);
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

