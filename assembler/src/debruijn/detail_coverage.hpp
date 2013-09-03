#pragma once

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

// PREVIOUS VERSION
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

		//FILE* file;
		double CountInCoverage ( const Sequence& seq, size_t size_bound ) const {

			unsigned len = (unsigned) seq.size();
			double edge_coverage_in(0.0);
			for (unsigned i = 0; i < size_bound; ++i) {
				runtime_k::RtSeq kmer_in(K_);
				for ( unsigned j = i; j < K_ + i && j < len; ++j) {
					kmer_in <<= seq[j];
				}
				edge_coverage_in += kmer_index_[kmer_in].count;
				//fprintf(file,"%d, ", kmer_index_[kmer_in].count);
			}
			//fprintf(file,"\n");
			return edge_coverage_in / (double) size_bound;
		}


		double CountOutCoverage ( const Sequence& seq, size_t size_bound ) const {

			unsigned len = (unsigned) seq.size();
			double edge_coverage_out(0.0);
			for (unsigned i = 0; i < size_bound; ++i) {
				runtime_k::RtSeq kmer_out(K_);
				for ( unsigned j = len-K_-i; j < len - i; ++j) {
					kmer_out <<= seq[j];
				}
				edge_coverage_out += kmer_index_[kmer_out].count;
			}
			return edge_coverage_out / (double) size_bound;
		}
	public:


		FlankingCoverage( const Graph& g, const Index& kmer_index, unsigned avg, unsigned K ) : kmer_index_(kmer_index), average_const_(avg), K_(K) {


			//file = fopen("/home/ksenia/coverage.log", "w");
  			for (auto e_iter = g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
				DEBUG( g.int_id(*e_iter));
				double edge_coverage_in(0.0), edge_coverage_out(0.0);
				//averageConst = graph.g.length(*e_iter);
				unsigned size_bound = g.length(*e_iter); //average_const_;
				if ( average_const_ > g.length(*e_iter) ){
					size_bound = (unsigned) g.length(*e_iter);
				}
				auto seq = g.EdgeNucls(*e_iter);
				edge_coverage_in = CountInCoverage( seq, size_bound );
				//fprintf(file,"%d: ", g.int_id(*e_iter));
				edge_coverage_out = CountOutCoverage( seq, size_bound );
				in_coverage_.insert(std::make_pair( *e_iter, edge_coverage_in));
				//inCoverage.insert(std::make_pair( *e_iter, graph.g.coverage(*e_iter)));
				out_coverage_.insert(std::make_pair( *e_iter, edge_coverage_out));
				//outCoverage.insert(std::make_pair( *e_iter, graph.g.coverage(*e_iter)));
			}
			//fclose(file);

		}
	private:
		DECL_LOGGER("FlankingCoverage");

	};
//
//template<class Graph, class Index>
//class FlankingCoverage : public GraphActionHandler<Graph> {
//  typedef GraphActionHandler<Graph> base;
//  typedef typename Graph::EdgeId EdgeId;
//  typedef typename Graph::VertexId VertexId;
//  /*
//   Iterates over kmer index saves values of coverage on the ends of edges
//   */
//
//  std::map<EdgeId, double> in_coverage_;
//  std::map<EdgeId, double> out_coverage_;
//  const Index& kmer_index_;
//  const size_t averaging_range_;
//
//  double CountAvgCoverage(EdgeId e, size_t offset) const {
//    size_t k = this->g().k();
//    VERIFY(offset == 0 || offset + averaging_range_ == this->g().length(e));
//    unsigned size_bound = std::min(averaging_range_, this->g().length(e));
//    const Sequence& seq = this->g().EdgeNucls(e);
//
//    size_t edge_coverage_in = 0;
//
//    runtime_k::RtSeq kpomer(k + 1, seq, offset);
//    kpomer >>= 0;
//    for (size_t i = 0; i < size_bound; ++i) {
//      kpomer <<= seq[offset + i + k];
//      VERIFY(kmer_index_.contains(kpomer));
//      edge_coverage_in += kmer_index_[kpomer].count;
//    }
//
//    return double(edge_coverage_in) / size_bound;
//  }
//
//  double CountInCoverage(EdgeId e) const {
//    return CountAvgCoverage(e, 0);
//  }
//
//  double CountOutCoverage(EdgeId e) const {
//    return CountAvgCoverage(e, std::max(this->g().length(e), averaging_range_) - averaging_range_);
//  }
//
// public:
//
//  FlankingCoverage(const Graph& g, const Index& kmer_index,
//                   unsigned averaging_range)
//      : base(g, "FlankingCoverage"),
//        kmer_index_(kmer_index),
//        averaging_range_(averaging_range) {
//    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
//      EdgeId e = *it;
//      in_coverage_.insert(std::make_pair(*it, CountInCoverage(e)));
//      out_coverage_.insert(std::make_pair(*it, CountOutCoverage(e)));
//    }
//  }
//
//  //todo rename
//  double GetInCov(EdgeId edge) const {
//    return get(in_coverage_, edge);
//  }
//
//  //todo rename
//  double GetOutCov(EdgeId edge) const {
//    return get(out_coverage_, edge);
//  }
//
//  /*virtual */
//  void HandleAdd(EdgeId e) {
//    in_coverage_.insert(std::make_pair(e, CountInCoverage(e)));
//    out_coverage_.insert(std::make_pair(e, CountOutCoverage(e)));
//  }
//
//  /*virtual*/
//  void HandleDelete(EdgeId e) {
//    in_coverage_.erase(e);
//    out_coverage_.erase(e);
//  }
//
//  double LocalCoverage(EdgeId e, VertexId v) const {
//    if (this->g().EdgeStart(e) == v) {
//      return GetInCov(e);
//    } else if (this->g().EdgeEnd(e) == v) {
//      return GetOutCov(e);
//    } else {
//      VERIFY(false);
//      return 0.0;
//    }
//  }
//
// private:
//  DECL_LOGGER("FlankingCoverage");
//
//};

}
