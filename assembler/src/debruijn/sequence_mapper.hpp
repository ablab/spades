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

    virtual size_t KmerSize() const = 0;
};

template<class Graph, class Index>
class NewExtendedSequenceMapper: public SequenceMapper<Graph> {

 using SequenceMapper<Graph>::g_;

 public:
  typedef std::vector<MappingRange> RangeMappings;
  typedef MappingPathFixer<Graph> GraphMappingPathFixer;

 private:
  const Index& index_;
  typedef typename Index::KMer Kmer;
  typedef KmerMapper<Graph, Kmer> KmerSubs;
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
