#ifndef DETAIL_COVERAGE_HPP
#define DETAIL_COVERAGE_HPP

#include "debruijn_kmer_index.hpp"
#include "graph_pack.hpp"
#include "verify.hpp"
#include "graphio.hpp"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>

namespace debruijn_graph {

	template <class IdType>
	struct FlankingCoverage{
	/*
	Iterates over kmer index saves values of coverage on the ends of edges
	*/

		std::map<EdgeId, double> inCoverage;
		std::map<EdgeId, double> outCoverage;
		int averageConst;
		
		FlankingCoverage( const conj_graph_pack& graph, const int avg ) : averageConst(avg) {

	
			size_t K = cfg::get().K + 1;
			DeBruijnEdgeIndex<EdgeId, runtime_k::RtSeq> kmerIndex(graph.index.inner_index().K(), cfg::get().output_dir);
			if (cfg::get().developer_mode) {

				std::string path;
				if (cfg::get().entry_point < ws_repeats_resolving) 
					path = cfg::get().output_dir + "/saves/debruijn_kmer_index_after_construction";
				else
					path = cfg::get().load_from + "/debruijn_kmer_index_after_construction";
				bool val = LoadEdgeIndex(path, kmerIndex);
				VERIFY_MSG(val, "can not open file "+path+".kmidx");
				INFO("Updating index from graph started");
				DeBruijnEdgeIndexBuilder<runtime_k::RtSeq>().UpdateIndexFromGraph(kmerIndex, graph.g);

			}
			else {

				//INFO("Updating index from graph started");
				//DeBruijnEdgeIndexBuilder<runtime_k::RtSeq>().UpdateIndexFromGraph(graph.index.inner_index(), graph.g);
				
			}

			int counter = 0;
  			for (auto e_iter = graph.g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
				
				auto seq = graph.g.EdgeNucls(*e_iter);
				int len = seq.size();
				double coverage_in(0.0), coverage_out(0.0);
				//averageConst = graph.g.length(*e_iter);
				int sizeBound = averageConst;
				
				if ( averageConst > graph.g.length(*e_iter) ){
					sizeBound = graph.g.length(*e_iter);		
				}

				// count incoming coverage
				for (int i = 0; i < sizeBound; ++i) {
					runtime_k::RtSeq kmer_in(K);
					for ( int j = i; j < K +i && j < len; ++j) {
						kmer_in <<= seq[j];
					}
					if (cfg::get().developer_mode) {
						if (kmerIndex.contains(kmer_in)) {
							coverage_in += kmerIndex[kmer_in].count_;
						}
					}
					else {
						if (graph.index.inner_index().contains(kmer_in)) {
							coverage_in += graph.index.inner_index()[kmer_in].count_;
						}
					}
				}
					
				coverage_in = coverage_in / sizeBound;

				// count outgoing coverage
				for (int i = 0; i < sizeBound; ++i) {
					runtime_k::RtSeq kmer_out;
					for ( int j = len-K-i; j >= 0 && j < len - i; ++j) {
							
						kmer_out <<= seq[j];
					}
					
					if (cfg::get().developer_mode) {
						if (kmerIndex.contains(kmer_out)) {
							coverage_out += kmerIndex[kmer_out].count_;
						}
					}
					else {
						if (graph.index.inner_index().contains(kmer_out)) {
							coverage_out += graph.index.inner_index()[kmer_out].count_;
						}
					}


				}
				coverage_out = coverage_out / sizeBound;
		
				//if (coverage_in != graph.g.coverage(*e_iter))
				//	std::cout << graph.g.int_id(*e_iter) << ": " << coverage_in << " " << coverage_out << " " << graph.g.coverage(*e_iter) << " " << graph.g.length(*e_iter) << std::endl;
				inCoverage.insert(std::make_pair( *e_iter, coverage_in));
				//inCoverage.insert(std::make_pair( *e_iter, graph.g.coverage(*e_iter)));
				outCoverage.insert(std::make_pair( *e_iter, coverage_out));
				//outCoverage.insert(std::make_pair( *e_iter, graph.g.coverage(*e_iter)));
				//	std::cout << graph.g.int_id(*e_iter) << " " << " " << graph.g.int_id( graph.g.EdgeStart(*e_iter)) << " " << graph.g.int_id(graph.g.EdgeEnd(*e_iter)) << " " << coverage_out << " " << graph.g.coverage(*e_iter) << " " << graph.g.length(*e_iter) << std::endl;
				//}
			}
			
			if (cfg::get().developer_mode) {
				kmerIndex.clear();
			}
		}
	};

}

#endif
