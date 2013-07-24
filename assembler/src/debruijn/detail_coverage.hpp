#ifndef DETAIL_COVERAGE_HPP
#define DETAIL_COVERAGE_HPP

#include "indices/debruijn_kmer_index.hpp"
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


	template <class Graph, class Index>
	class FlankingCoverage{
		typedef typename Graph::EdgeId EdgeId;
		typedef typename Graph::VertexId VertexId;
//		typedef DeBruijnEdgeIndex<KmerStoringDeBruijnEdgeIndex<Graph, runtime_k::RtSeq>> Index;
	/*
	Iterates over kmer index saves values of coverage on the ends of edges
	*/

		std::map<EdgeId, double> in_coverage_;
		std::map<EdgeId, double> out_coverage_;
		const Index& kmer_index_;
		const unsigned average_const_;
		const unsigned K_;

	public:

		double GetInCov(EdgeId edge ) const {
			auto coverage = in_coverage_.find(edge);
			VERIFY(coverage != in_coverage_.end());
			return coverage->second;
		}

		double GetOutCov(EdgeId edge) const {

			auto coverage = out_coverage_.find(edge);
			VERIFY(coverage != out_coverage_.end());
			return coverage->second;
		}

	private:

		double CountInCoverage ( const Sequence& seq, size_t size_bound ) const {

			unsigned len = (unsigned) seq.size();
			double edge_coverage_in(0.0);

			for (unsigned i = 0; i < size_bound; ++i) {
				runtime_k::RtSeq kmer_in(K_);
				for ( unsigned j = i; j < K_ + i && j < len; ++j) {
					kmer_in <<= seq[j];
				}

				edge_coverage_in += kmer_index_[kmer_in].count;
			}

			return edge_coverage_in / (double) size_bound;


		}

		double CountOutCoverage ( const Sequence& seq, size_t size_bound ) const {

			unsigned len = (unsigned) seq.size();
			double edge_coverage_out(0.0);

			for (unsigned i = 0; i < size_bound; ++i) {
				runtime_k::RtSeq kmer_out;
				for ( unsigned j = len-K_-i; j < len - i; ++j) {
							
					kmer_out <<= seq[j];
				}

				edge_coverage_out += kmer_index_[kmer_out].count;


				}
			return edge_coverage_out / (double) size_bound;


		}
	public:


		FlankingCoverage( const Graph& g, const Index& kmer_index, unsigned avg, unsigned K ) : kmer_index_(kmer_index), average_const_(avg), K_(K) {


  			for (auto e_iter = g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {

				double edge_coverage_in(0.0), edge_coverage_out(0.0);
				//averageConst = graph.g.length(*e_iter);
				unsigned size_bound = average_const_;

				if ( average_const_ > g.length(*e_iter) ){
					size_bound = (unsigned) g.length(*e_iter);
				}

				auto seq = g.EdgeNucls(*e_iter);
				edge_coverage_in = CountInCoverage( seq, size_bound );

				edge_coverage_out = CountOutCoverage( seq, size_bound );
				in_coverage_.insert(std::make_pair( *e_iter, edge_coverage_in));
				//inCoverage.insert(std::make_pair( *e_iter, graph.g.coverage(*e_iter)));
				out_coverage_.insert(std::make_pair( *e_iter, edge_coverage_out));
				//outCoverage.insert(std::make_pair( *e_iter, graph.g.coverage(*e_iter)));
			}

		}
	private:
		DECL_LOGGER("FlankingCoverage");

	};

}

#endif
