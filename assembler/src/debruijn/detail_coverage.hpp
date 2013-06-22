//#ifndef DETAIL_COVERAGE_HPP
//#define DETAIL_COVERAGE_HPP
//
//#include "indices/debruijn_kmer_index.hpp"
//#include "graph_pack.hpp"
//#include "verify.hpp"
//#include "graphio.hpp"
//#include <vector>
//#include <map>
//#include <set>
//#include <string>
//#include <iostream>
//#include <fstream>
//
//namespace debruijn_graph {
//
//	template <class IdType>
//	struct FlankingCoverage{
//	/*
//	Iterates over kmer index saves values of coverage on the ends of edges
//	*/
//
//		std::map<EdgeId, int> inCoverage;
//		std::map<EdgeId, int> outCoverage;
//		static int averageConst;
//
//		FlankingCoverage( const conj_graph_pack& graph ) {
//
//			size_t K = cfg::get().K + 1;
//			DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> kmerIndex(graph.index.inner_index().K(), cfg::get().output_dir);
//			std::string path = cfg::get().output_dir + "/saves/debruijn_kmer_index_after_construction";
//			bool val = LoadEdgeIndex(path, kmerIndex);
//			VERIFY_MSG(val, "can not open file "+path+".kmidx");
//
//
//			INFO("Updating index from graph started");
//			DeBruijnEdgeIndexBuilder<runtime_k::RtSeq>().UpdateIndexFromGraph(kmerIndex, graph.g);
//
//			int counter = 0;
//  			for (auto e_iter = graph.g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
//
//				auto seq = graph.g.EdgeNucls(*e_iter);
//				int len = seq.size();
//				double coverage_in(0.0), coverage_out(0.0);
//				//averageConst = graph.g.length(*e_iter);
//				int currentNumberOfIterations = averageConst;
//
//				if ( averageConst > graph.g.length(*e_iter) ){
//					currentNumberOfIterations = graph.g.length(*e_iter);
//				}
//
//				// count incoming coverage
//				for (int i = 0; i < currentNumberOfIterations; ++i) {
//					runtime_k::RtSeq kmer_in(K);
//					for ( int j = i; j < K +i && j < len; ++j) {
//						kmer_in <<= seq[j];
//					}
//					if (!graph.index.inner_index().contains(kmer_in)) {
//						continue;
//					}
//					coverage_in += graph.index.inner_index()[kmer_in].count_;
//				}
//
//				coverage_in = coverage_in / currentNumberOfIterations;
//
//				// count outgoing coverage
//				for (int i = 0; i < currentNumberOfIterations; ++i) {
//					runtime_k::RtSeq kmer_out;
//					for ( int j = len-K-i; j >= 0 && j < len - i; ++j) {
//
//						kmer_out <<= seq[j];
//					}
//					if (!graph.index.inner_index().contains(kmer_out)){
//						continue;
//					}
//					coverage_out += graph.index.inner_index()[kmer_out].count_;
//
//				}
//				coverage_out = coverage_out / currentNumberOfIterations;
//
//				inCoverage.insert(std::make_pair( *e_iter, coverage_in));
//				//inCoverage.insert(std::make_pair( *e_iter, graph.g.coverage(*e_iter)));
//				outCoverage.insert(std::make_pair( *e_iter, coverage_out));
//				//outCoverage.insert(std::make_pair( *e_iter, graph.g.coverage(*e_iter)));
//				if (coverage_out - graph.g.coverage(*e_iter) > 1) {
//					counter++;
//					std::cout << graph.g.int_id(*e_iter) << " " << coverage_out << " " << graph.g.coverage(*e_iter) << std::endl;
//				}
//			}
//
//			std::cout << "COUNTER: " << counter << std::endl;
//
//		}
//	};
//
//	template <class IdType>
//	int FlankingCoverage<IdType>::averageConst = 50;
//}
//
//#endif
