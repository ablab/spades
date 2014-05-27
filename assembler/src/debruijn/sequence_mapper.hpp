//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "omni/omni_utils.hpp"
#include "sequence/sequence_tools.hpp"
#include "omni/path_processor.hpp"

#include "runtime_k.hpp"
#include "edge_index.hpp"
#include "kmer_mapper.hpp"

#include <cstdlib>

namespace debruijn_graph {
///**
// * This class finds how certain sequence is mapped to genome. As it is now it works correct only if sequence
// * is mapped to graph ideally and in unique way.
// */
//template <class Graph, class Index>
//class SimpleSequenceMapper {
//  /*
//};
//
//template<class Graph>
//class SimpleSequenceMapper<Graph, runtime_k::RtSeq> {
//*/
// public:
//  typedef typename Graph::EdgeId EdgeId;
//  typedef typename Index::KMer Kmer;
//
// private:
//  const Graph& g_;
//  const Index &index_;
//  size_t k_;
//
//  bool TryThread(Kmer &kmer, std::vector<EdgeId> &passed,
//                 size_t& endPosition) const {
//    VERIFY(passed.size() > 0);
//    EdgeId last = passed[passed.size() - 1];
//    if (endPosition + 1 < g_.length(last)) {
//      if (g_.EdgeNucls(last)[endPosition + k_] == kmer[k_ - 1]) {
//        endPosition++;
//        return true;
//      }
//    } else {
//      VertexId v = g_.EdgeEnd(last);
//      for (auto I = g_.out_begin(v), E = g_.out_end(v); I != E; ++I) {
//        EdgeId edge = *I;
//        if (g_.EdgeNucls(edge)[k_ - 1] == kmer[k_ - 1]) {
//          passed.push_back(edge);
//          endPosition = 0;
//          return true;
//        }
//      }
//    }
//    return false;
//  }
//
//  bool FindKmer(Kmer& kmer, std::vector<EdgeId> &passed, size_t& startPosition,
//                size_t& endPosition) const {
//    //TRACE("CONTAINS kmer " << " " << omp_get_thread_num() );
//    pair<EdgeId, size_t> position = index_.get(kmer);
//    if (position.second != -1u) {
//      //TRACE("YES CONTAINS " << omp_get_thread_num());
//      //DEBUG("LENGTH " << g_.length(position.first));
//      endPosition = position.second;
//      if (passed.empty()) {
//        startPosition = position.second;
//      }
//      if (passed.empty() || passed.back() != position.first) {
//        passed.push_back(position.first);
//      }
//
//      return true;
//    }
//    return false;
//  }
//
//  bool ProcessKmer(Kmer &kmer, std::vector<EdgeId> &passed, size_t &startPosition,
//                   size_t &endPosition, bool valid) const {
//    //DEBUG("process kmer started " << omp_get_thread_num() << " valid " << valid);
//    if (valid) {
//      return TryThread(kmer, passed, endPosition);
//    } else {
//      return FindKmer(kmer, passed, startPosition, endPosition);
//    }
//    return false;
//    //DEBUG("process kmer finished " << omp_get_thread_num());
//  }
//
// public:
//  /**
//   * Creates SimpleSequenceMapper for given graph. Also requires index_ which should be synchronized
//   * with graph.
//   * @param g graph sequences should be mapped to
//   * @param index index synchronized with graph
//   */
//  SimpleSequenceMapper(const Graph& g, const Index& index) :
//      g_(g), index_(index), k_(g.k()+1) {
//  }
//
//  /**
//   * Finds a path in graph which corresponds to given sequence.
//   * @read sequence to be mapped
//   */
//
//  Path<EdgeId> MapSequence(const Sequence &read) const {
//    std::vector<EdgeId> passed;
//    //TRACE("Mapping sequence");
//    if (read.size() <= k_ - 1) {
//      return Path<EdgeId>();
//    }
//
//    Kmer kmer = read.start<Kmer>(k_);
//    //DEBUG("started " << kmer.str() << omp_get_thread_num() );
//    size_t startPosition = -1ul;
//    size_t endPosition = -1ul;
//    bool valid = ProcessKmer(kmer, passed, startPosition, endPosition,
//                             false);
//    for (size_t i = k_; i < read.size(); ++i) {
//      kmer <<= read[i];
//      //DEBUG("shifted " << kmer.str() << omp_get_thread_num());
//      valid = ProcessKmer(kmer, passed, startPosition, endPosition,
//                          valid);
//    }
//    //DEBUG("Path got " << omp_get_thread_num());
//    Path<EdgeId> ans(passed, startPosition, endPosition + 1);
//    return ans;
//  }
//
//};

//template<class Graph, class Seq = runtime_k::RtSeq>
//class ExtendedSequenceMapper {
// public:
//  typedef typename Graph::EdgeId EdgeId;
//  typedef std::vector<MappingRange> RangeMappings;
//  typedef Seq Kmer;
//  typedef EdgeIndex<Graph, Kmer> Index;
//  typedef KmerMapper<Graph, Kmer> KmerSubs;
//
// private:
//  const Graph& g_;
//  const Index& index_;
//  const KmerSubs& kmer_mapper_;
//  size_t k_;
//
//  void FindKmer(Kmer kmer, size_t kmer_pos, std::vector<EdgeId> &passed,
//                RangeMappings& range_mappings) const {
//
//    if (index_.contains(kmer)) {
//      pair<EdgeId, size_t> position = index_.get(kmer);
//      if (passed.empty() || passed.back() != position.first
//          || kmer_pos != range_mappings.back().initial_range.end_pos
//          || position.second + 1
//          < range_mappings.back().mapped_range.end_pos) {
//        passed.push_back(position.first);
//        MappingRange mapping_range(Range(kmer_pos, kmer_pos + 1),
//                                   Range(position.second, position.second + 1));
//        range_mappings.push_back(mapping_range);
//      } else {
//        range_mappings.back().initial_range.end_pos = kmer_pos + 1;
//        range_mappings.back().mapped_range.end_pos = position.second
//                                                     + 1;
//      }
//    }
//  }
//
//  void ProcessKmer(Kmer kmer, size_t kmer_pos, std::vector<EdgeId> &passed,
//                   RangeMappings& interval_mapping) const {
//    kmer = kmer_mapper_.Substitute(kmer);
//    FindKmer(kmer, kmer_pos, passed, interval_mapping);
//  }
//
// public:
//  ExtendedSequenceMapper(const Graph& g,
//                         const Index& index,
//                         const KmerSubs& kmer_mapper,
//                         size_t k) :
//      g_(g), index_(index), kmer_mapper_(kmer_mapper), k_(k) {
//  }
//
//  MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
//    std::vector<EdgeId> passed_edges;
//    RangeMappings range_mapping;
//
//    if (sequence.size() < k_) {
//      return MappingPath<EdgeId>();
//    }
//    Kmer kmer = sequence.start<Kmer>(k_);
//    //kmer >>= 0;
//    ProcessKmer(kmer, 0, passed_edges, range_mapping);
//    for (size_t i = k_; i < sequence.size(); ++i) {
//      kmer <<= sequence[i];
//      ProcessKmer(kmer, i - k_ + 1, passed_edges, range_mapping);
//    }
//
//    //DEBUG
//    //		for (size_t i = 0; i < passed_edges.size(); ++i) {
//    //			cerr << passed_edges[i].int_id() << " (" << range_mapping[i] << ")"<< "; ";
//    //		}
//    //		cerr << endl;
//    //DEBUG
//
//    return MappingPath<EdgeId>(passed_edges, range_mapping);
//  }
//};



template<class Graph>
class SequenceMapper {
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef runtime_k::RtSeq Kmer;

protected:
    const Graph& g_;

public:
    SequenceMapper(const Graph& g): g_(g) {

    }

    virtual ~SequenceMapper() {

    }

    virtual MappingPath<EdgeId> MapSequence(const Sequence &sequence) const = 0;

  
    MappingPath<EdgeId> MapRead(const io::SingleRead &read) const {
      VERIFY(read.IsValid());
      return MapSequence(read.sequence());
    }

    virtual pair<EdgeId, size_t> GetKmerPos(const Kmer& kmer) const = 0;

    virtual size_t KmerSize() const = 0;
};

//todo compare performance
template<class Graph, class Index>
class NewExtendedSequenceMapper: public SequenceMapper<Graph> {

 using SequenceMapper<Graph>::g_;

 public:
  typedef typename SequenceMapper<Graph>::Kmer Kmer;
  typedef std::vector<MappingRange> RangeMappings;
  typedef KmerMapper<Graph, Kmer> KmerSubs;
  typedef MappingPathFixer<Graph> GraphMappingPathFixer;

 private:
  const Index& index_;
  const KmerSubs& kmer_mapper_;
  const GraphMappingPathFixer path_fixer_;
  size_t k_;
  bool optimization_on_;
  //	mutable size_t mapped_;
  //	mutable size_t unmapped_;

  bool FindKmer(const Kmer &kmer, size_t kmer_pos, std::vector<EdgeId> &passed,
                RangeMappings& range_mappings) const {
    std::pair<EdgeId, size_t> position = index_.get(kmer);
    if (position.second != -1u/*index contains this k-mer*/) {
      if (passed.empty() || passed.back() != position.first ||
          kmer_pos != range_mappings.back().initial_range.end_pos ||
          position.second + 1 < range_mappings.back().mapped_range.end_pos) {
        passed.push_back(position.first);
        range_mappings.push_back(
            MappingRange(Range(kmer_pos, kmer_pos + 1),
                         Range(position.second, position.second + 1)));
      } else {
        range_mappings.back().initial_range.end_pos = kmer_pos + 1;
        range_mappings.back().mapped_range.end_pos = position.second + 1;
      }
      return true;
    }
    return false;
  }

  bool TryThread(const Kmer& kmer, size_t kmer_pos, std::vector<EdgeId> &passed,
                 RangeMappings& range_mappings) const {
    EdgeId last_edge = passed.back();
    size_t end_pos = range_mappings.back().mapped_range.end_pos;
    if (end_pos < g_.length(last_edge)) {
      if (g_.EdgeNucls(last_edge)[end_pos + k_ - 1] == kmer[k_ - 1]) {
        range_mappings.back().initial_range.end_pos++;
        range_mappings.back().mapped_range.end_pos++;
        return true;
      }
    } else {
      VertexId v = g_.EdgeEnd(last_edge);

      if(!optimization_on_)
    	  if(g_.OutgoingEdgeCount(v) > 1)
    		  return false;

      for (auto I = g_.out_begin(v), E = g_.out_end(v); I != E; ++I) {
        EdgeId edge = *I;
        if (g_.EdgeNucls(edge)[k_ - 1] == kmer[k_ - 1]) {
          passed.push_back(edge);
          range_mappings.push_back(
              MappingRange(Range(kmer_pos, kmer_pos + 1),
                           Range(0, 1)));
          return true;
        }
      }
    }
    return false;
  }

  bool Substitute(Kmer& kmer) const {
    Kmer subs = kmer_mapper_.Substitute(kmer);
    if (subs != kmer) {
      kmer = subs;
      return true;
    }
    return false;
  }

  bool ProcessKmer(Kmer kmer, size_t kmer_pos, std::vector<EdgeId> &passed_edges,
                   RangeMappings& range_mapping, bool try_thread) const {
    if (try_thread) {
      if (!TryThread(kmer, kmer_pos, passed_edges, range_mapping)) {
        Substitute(kmer);
        FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
        return false;
      } else {
        return true;
      }
    } else {
      if (!Substitute(kmer)) {
        return FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
      } else {
        FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
        return false;
      }
    }
    //		if (!Substitute(kmer)) {
    //			if (try_thread) {
    //				return TryThread(kmer, kmer_pos, passed_edges, range_mapping);
    //			} else {
    //				return FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
    //			}
    //		} else {
    //			FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
    //			return false;
    //		}
  }

 public:
  NewExtendedSequenceMapper(const Graph& g,
                            const Index& index,
                            const KmerSubs& kmer_mapper,
			    bool optimization_on = true) :
      SequenceMapper<Graph>(g), index_(index), kmer_mapper_(kmer_mapper), path_fixer_(g), k_(g.k()+1),
	optimization_on_(optimization_on) { }

  ~NewExtendedSequenceMapper() {
    //		TRACE("In destructor of sequence mapper");
    //		TRACE(mapped_ << " sequences were mapped");
    //		TRACE(unmapped_ << " sequences couldn't be mapped");
  }

  MappingPath<EdgeId> MapSequence(const Sequence &sequence) const {
    std::vector<EdgeId> passed_edges;
    RangeMappings range_mapping;

    if (sequence.size() < k_) {
      return MappingPath<EdgeId>();
    }

    Kmer kmer = sequence.start<Kmer>(k_);
    //kmer >>= 0;
    bool try_thread = false;
    try_thread = ProcessKmer(kmer, 0, passed_edges,
                             range_mapping, try_thread);
    for (size_t i = k_; i < sequence.size(); ++i) {
      kmer <<= sequence[i];
      try_thread = ProcessKmer(kmer, i - k_ + 1, passed_edges,
                               range_mapping, try_thread);
    }

    //		if (passed_edges.empty()) {
    ////			TRACE("Sequence " << sequence << "couldn't be mapped");
    //			unmapped_++;
    //			//todo maybe check path consistency?
    //		} else {
    //			mapped_++;
    //		}

    return MappingPath<EdgeId>(passed_edges, range_mapping);
  }

  pair<EdgeId, size_t> GetKmerPos(const Kmer& kmer) const {
    VERIFY(kmer.size() == k_);

    return index_.get(kmer_mapper_.Substitute(kmer));
  }

  size_t KmerSize() const {
      return k_;
  }

  vector<EdgeId> FindReadPath(const MappingPath<EdgeId>& mapping_path) const {
        if (!IsMappingPathValid(mapping_path)) {
            TRACE("read unmapped");
            return vector<EdgeId>();
        }
        vector<EdgeId> corrected_path = path_fixer_.DeleteSameEdges(
                mapping_path.simple_path());
        vector<EdgeId> fixed_path = path_fixer_.TryFixPath(corrected_path);
        if (!path_fixer_.CheckContiguous(fixed_path)) {
            TRACE("read unmapped");
            std::stringstream debug_stream;
            for (size_t i = 0; i < fixed_path.size(); ++i) {
                debug_stream << g_.int_id(fixed_path[i]) << " ";
            }TRACE(debug_stream.str());
            return vector<EdgeId>();
        }
        return fixed_path;
    }

private:
    bool IsMappingPathValid(const MappingPath<EdgeId>& path) const {
        return path.size() != 0;
    }
    DECL_LOGGER("NewExtendedSequenceMapper");
};


template<class gp_t>
std::shared_ptr<NewExtendedSequenceMapper<typename gp_t::graph_t, typename gp_t::index_t> > MapperInstance(const gp_t& gp) {
  return std::make_shared<NewExtendedSequenceMapper<typename gp_t::graph_t, typename gp_t::index_t> >(gp.g, gp.index, gp.kmer_mapper);
}


}
