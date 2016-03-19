//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "sequence/sequence_tools.hpp"
#include "simplification/path_processor.hpp"
#include "graph_support/basic_graph_stats.hpp"

#include "sequence/runtime_k.hpp"
#include "edge_index.hpp"
#include "kmer_mapper.hpp"

#include <cstdlib>
#include "graph_support/basic_graph_stats.hpp"

namespace debruijn_graph {
using omnigraph::MappingPath;
using omnigraph::Path;
using omnigraph::MappingRange;
using omnigraph::Range;

template<class Graph>
MappingPath<typename Graph::EdgeId> ConjugateMapping(const Graph& g, 
                                              const MappingPath<typename Graph::EdgeId>& mp, 
                                              size_t sequence_length) {
    MappingPath<typename Graph::EdgeId> answer;
    for (size_t i = mp.size(); i > 0; --i) {
        auto p = mp[i-1];
        auto e = p.first;
        MappingRange mr = p.second;
        answer.push_back(g.conjugate(e), 
                        MappingRange(mr.initial_range.Invert(sequence_length - g.k()), 
                        mr.mapped_range.Invert(g.length(e))));
    }
    return answer;
}

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
//      VERIFY(read.IsValid());
        DEBUG(read.name() << " is mapping");
        string s = read.GetSequenceString();
        size_t l = 0, r = 0;
        MappingPath<EdgeId> result;
        for(size_t i = 0; i < s.size(); i++) {
            if (read.GetSequenceString()[i] == 'N') {
                if (r > l) {
                    result.join(MapSequence(Sequence(s.substr(l, r - l))), int(l));
                }
                r = i + 1;
                l = i + 1;
            } else {
                r++;
            }
        }
        if (r > l) {
            result.join(MapSequence(Sequence(s.substr(l, r - l))), int(l));
        }
        DEBUG(read.name() << " is mapped");
        DEBUG("Number of edges is " << result.size());

      return result;
    }

    virtual size_t KmerSize() const = 0;
};

template<class Graph>
class MappingPathFixer {
public:

    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    MappingPathFixer(const Graph& graph)
            : g_(graph) {
    }

    bool CheckContiguous(const vector<typename Graph::EdgeId>& path) const {
        for (size_t i = 1; i < path.size(); ++i) {
            if (g_.EdgeEnd(path[i - 1]) != g_.EdgeStart(path[i]))
                return false;
        }
        return true;
    }

    Path<EdgeId> TryFixPath(const Path<EdgeId>& path, size_t length_bound = 70) const {
        return Path<EdgeId>(TryFixPath(path.sequence(), length_bound), path.start_pos(), path.end_pos());
    }

    vector<EdgeId> TryFixPath(const vector<EdgeId>& edges, size_t length_bound = 70) const {
        vector<EdgeId> answer;
        if (edges.empty()) {
            //          WARN("Mapping path was empty");
            return vector<EdgeId>();
        }
        answer.push_back(edges[0]);
        for (size_t i = 1; i < edges.size(); ++i) {
            if (g_.EdgeEnd(edges[i - 1]) != g_.EdgeStart(edges[i])) {
                vector<EdgeId> closure = TryCloseGap(g_.EdgeEnd(edges[i - 1]),
                                                     g_.EdgeStart(edges[i]),
                                                     length_bound);
                answer.insert(answer.end(), closure.begin(), closure.end());
            }
            answer.push_back(edges[i]);
        }
        return answer;
    }

    vector<EdgeId> DeleteSameEdges(const vector<EdgeId>& path) const {
        vector<EdgeId> result;
        if (path.empty()) {
            return result;
        }
        result.push_back(path[0]);
        for (size_t i = 1; i < path.size(); ++i) {
            if (path[i] != result[result.size() - 1]) {
                result.push_back(path[i]);
            }
        }
        return result;
    }

private:
    vector<EdgeId> TryCloseGap(VertexId v1, VertexId v2, size_t length_bound) const {
        if (v1 == v2)
            return vector<EdgeId>();
        TRACE("Trying to close gap between v1=" << g_.int_id(v1) << " and v2=" << g_.int_id(v2));
        omnigraph::PathStorageCallback<Graph> path_store(g_);

        TRACE("Path storage callback created");
        //todo reduce value after investigation
        omnigraph::ProcessPaths(g_, 0, length_bound, v1, v2, path_store);

        TRACE("Paths processed");
        if (path_store.size() == 0) {
            TRACE("Failed to find closing path");
            //          TRACE("Failed to close gap between v1=" << graph_.int_id(v1)
            //                          << " (conjugate "
            //                          << graph_.int_id(g_.conjugate(v1))
            //                          << ") and v2=" << g_.int_id(v2)
            //                          << " (conjugate "
            //                          << g_.int_id(g_.conjugate(v2)) << ")");
            //          return boost::none;
            return vector<EdgeId>();
        } else if (path_store.size() == 1) {
            TRACE("Unique closing path found");
        } else {
            TRACE("Several closing paths found, first chosen");
        }
        TRACE("Taking answer    ");
        vector<EdgeId> answer = path_store.paths().front();
        TRACE("Gap closed");
        TRACE( "Cumulative closure length is " << CumulativeLength(g_, answer));
        return answer;
    }
    const Graph& g_;
};

template<class Graph>
class ReadPathFinder {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph& g_;
    typedef MappingPathFixer<Graph> GraphMappingPathFixer;
    const GraphMappingPathFixer path_fixer_;
public:
    ReadPathFinder (const Graph& g) :
        g_(g), path_fixer_(g)
    {   }

    vector<EdgeId> FindReadPath(const MappingPath<EdgeId>& mapping_path) const {
          if (!IsMappingPathValid(mapping_path)) {
              TRACE("read unmapped");
              return vector<EdgeId>();
          }
          vector<EdgeId> corrected_path = path_fixer_.DeleteSameEdges(
                  mapping_path.simple_path());
          PrintPathInfo(corrected_path);
          if(corrected_path.size() != mapping_path.simple_path().size()) {
              DEBUG("Some edges were deleted");
          }
          vector<EdgeId> fixed_path = path_fixer_.TryFixPath(corrected_path);
          if (!path_fixer_.CheckContiguous(fixed_path)) {
              TRACE("read unmapped");
              std::stringstream debug_stream;
              for (size_t i = 0; i < fixed_path.size(); ++i) {
                  debug_stream << g_.int_id(fixed_path[i]) << " ";
              }
              TRACE(debug_stream.str());
              return vector<EdgeId>();
          } else {
              DEBUG("Path fix works");
          }
          return fixed_path;
      }

    vector<vector<EdgeId>> FindReadPathWithGaps(const MappingPath<EdgeId>& mapping_path) const {
          if (!IsMappingPathValid(mapping_path)) {
              TRACE("read unmapped");
              return vector<vector<EdgeId>>();
          }
          vector<EdgeId> corrected_path = path_fixer_.DeleteSameEdges(
                  mapping_path.simple_path());
          PrintPathInfo(corrected_path);
          if(corrected_path.size() != mapping_path.simple_path().size()) {
              DEBUG("Some edges were deleted");
          }
          vector<EdgeId> fixed_path = path_fixer_.TryFixPath(corrected_path);
          return SplitUnfixedPoints(fixed_path);
      }

private:

      vector<vector<EdgeId>> SplitUnfixedPoints(vector<EdgeId>& path) const {
          vector<vector<EdgeId>> result;
          size_t prev_start = 0;
          for (size_t i = 1; i < path.size(); ++i) {
              if (g_.EdgeEnd(path[i - 1]) != g_.EdgeStart(path[i])) {
                      result.push_back(vector<EdgeId>(path.begin() + prev_start, path.begin() + i));
                      prev_start = i;
              }
          }
          result.push_back(vector<EdgeId>(path.begin() + prev_start, path.end()));
          return result;
      }

      bool IsTip(VertexId v) const {
          return g_.IncomingEdgeCount(v) + g_.OutgoingEdgeCount(v) == 1;
      }

      bool IsMappingPathValid(const MappingPath<EdgeId>& path) const {
          return path.size() != 0;
      }

      void PrintPathInfo(vector<EdgeId>& corrected_path) const {
          for(size_t i = 0; i < corrected_path.size(); ++i) {
              DEBUG(i + 1 << "-th edge is " << corrected_path[i].int_id());
          }
      }
};

template<class Graph, class Index>
class NewExtendedSequenceMapper: public SequenceMapper<Graph> {

 using SequenceMapper<Graph>::g_;

 public:
  typedef std::vector<MappingRange> RangeMappings;

 private:
  const Index& index_;
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  typedef typename Index::KMer Kmer;
  typedef KmerMapper<Graph, Kmer> KmerSubs;
  const KmerSubs& kmer_mapper_;
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
      SequenceMapper<Graph>(g), index_(index), kmer_mapper_(kmer_mapper), k_(g.k()+1),
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

  DECL_LOGGER("NewExtendedSequenceMapper");
};


template<class gp_t>
std::shared_ptr<NewExtendedSequenceMapper<typename gp_t::graph_t, typename gp_t::index_t> > MapperInstance(const gp_t& gp) {
  return std::make_shared<NewExtendedSequenceMapper<typename gp_t::graph_t, typename gp_t::index_t> >(gp.g, gp.index, gp.kmer_mapper);
}

}
