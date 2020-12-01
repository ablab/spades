//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/basic_graph_stats.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/paths/path_processor.hpp"
#include "io/reads/single_read.hpp"

#include "sequence/sequence_tools.hpp"
#include "pipeline/graph_pack.hpp"

#include "kmer_mapper.hpp"
#include "edge_index.hpp"

#include <cstdlib>

namespace debruijn_graph {
using omnigraph::MappingPath;
using omnigraph::Path;
using omnigraph::MappingRange;

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
    typedef RtSeq Kmer;

    virtual ~SequenceMapper() = default;

    virtual MappingPath<EdgeId> MapSequence(const Sequence &sequence,
                                            bool only_simple = false) const = 0;

    virtual MappingPath<EdgeId> MapRead(const io::SingleRead &read,
                                        bool only_simple = false) const = 0;
};

template<class Graph>
class AbstractSequenceMapper : public SequenceMapper<Graph> {
protected:
    const Graph& g_;

public:
    AbstractSequenceMapper(const Graph& g)
            : g_(g) {}

    MappingPath<EdgeId> MapRead(const io::SingleRead &read,
                                bool only_simple = false) const override {
//      VERIFY(read.IsValid());
        DEBUG(read.name() << " is mapping");
        auto s = read.GetSequenceString();
        size_t l = 0, r = 0;
        MappingPath<EdgeId> result;
        for (size_t i = 0; i < s.size(); i++) {
            if (read.GetSequenceString()[i] == 'N') {
                if (r > l) {
                    result.join(this->MapSequence(Sequence(s.substr(l, r - l))), int(l));
                }
                r = i + 1;
                l = i + 1;
            } else {
                r++;
            }
        }
        if (r > l) {
            result.join(this->MapSequence(Sequence(s.substr(l, r - l))), int(l));
        }

        // FIXME: exit earlier
        if (only_simple && result.size() > 1)
            result = MappingPath<EdgeId>();
        
        DEBUG(read.name() << " is mapped, only simple mode: " << only_simple);
        DEBUG("Number of edges is " << result.size());

        return result;
    }
};

template<class Graph>
class DelegatingSequenceMapper : public SequenceMapper<Graph> {
public:
    typedef std::function<MappingPath<EdgeId> (const MappingPath<EdgeId>&, size_t)> ProcessingF;
private:
    std::shared_ptr<SequenceMapper<Graph>> inner_mapper_;
    ProcessingF processing_f_;

public:
    DelegatingSequenceMapper(std::shared_ptr<SequenceMapper<Graph>> inner_mapper,
                             ProcessingF processing_f) :
            inner_mapper_(inner_mapper), processing_f_(processing_f) {}

    MappingPath<EdgeId> MapSequence(const Sequence &s,
                                    bool only_simple = false) const override {
        return processing_f_(inner_mapper_->MapSequence(s, only_simple), s.size());
    }

    MappingPath<EdgeId> MapRead(const io::SingleRead &r,
                                bool only_simple = false) const override {
        return processing_f_(inner_mapper_->MapRead(r, only_simple), r.size());
    }
};

template<class Graph>
bool SpuriousMappingFilter(const Graph& /*g*/,
                           const MappingPath<EdgeId>& mapping_path,
                           size_t read_length,
                           size_t max_range,
                           size_t min_flank) {
    if (mapping_path.size() != 1)
        return true;
    
    Range read_range = mapping_path[0].second.initial_range;
    if (read_range.size() <= max_range &&
        read_range.start_pos >= min_flank &&
        read_range.end_pos + min_flank <= read_length)
        return false;

    return true;
}

template<class Graph>
bool CheckContiguous(const Graph &g, const std::vector<typename Graph::EdgeId> &path) {
    for (size_t i = 1; i < path.size(); ++i) {
        if (g.EdgeEnd(path[i - 1]) != g.EdgeStart(path[i]))
            return false;
    }
    return true;
}


template<class Graph>
class MappingPathFixer {
public:

    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    constexpr static size_t LENGTH_BOUND_DEFAULT = 70;

    MappingPathFixer(const Graph& graph)
            : g_(graph) {
    }

    Path<EdgeId> TryFixPath(const Path<EdgeId>& path, size_t length_bound = LENGTH_BOUND_DEFAULT) const {
        return Path<EdgeId>(TryFixPath(path.sequence(), length_bound), path.start_pos(), path.end_pos());
    }

    std::vector<EdgeId> TryFixPath(const std::vector<EdgeId>& edges, size_t length_bound = LENGTH_BOUND_DEFAULT) const {
        std::vector<EdgeId> answer;
        if (edges.empty()) {
            //          WARN("Mapping path was empty");
            return {};
        }
        answer.push_back(edges[0]);
        for (size_t i = 1; i < edges.size(); ++i) {
            if (g_.EdgeEnd(edges[i - 1]) != g_.EdgeStart(edges[i])) {
                auto closure = TryCloseGap(g_.EdgeEnd(edges[i - 1]), g_.EdgeStart(edges[i]), length_bound);
                answer.insert(answer.end(), closure.begin(), closure.end());
            }
            answer.push_back(edges[i]);
        }
        return answer;
    }

    std::vector<EdgeId> DeleteSameEdges(const std::vector<EdgeId>& path) const {
        std::vector<EdgeId> result;
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

    Path<EdgeId> DeleteSameEdges(const Path<EdgeId>& path) const {
        return Path<EdgeId>(DeleteSameEdges(path.sequence()), path.start_pos(), path.end_pos());
    }

    std::vector<EdgeId> TryCloseGap(VertexId v1, VertexId v2, size_t length_bound = LENGTH_BOUND_DEFAULT) const {
        if (v1 == v2)
            return {};
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
            return {};
        } else if (path_store.size() == 1) {
            TRACE("Unique closing path found");
        } else {
            TRACE("Several closing paths found, first chosen");
        }
        TRACE("Taking answer    ");
        const auto& answer = path_store.paths().front();
        TRACE("Gap closed");
        TRACE( "Cumulative closure length is " << CumulativeLength(g_, answer));
        return answer;
    }

private:
    const Graph& g_;
};

template<class Graph>
class ReadPathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& g_;
    const MappingPathFixer<Graph> path_fixer_;
    const bool skip_unfixed_;
    const size_t max_gap_fill_;

public:
    ReadPathFinder(const Graph& g, bool skip_unfixed = true, size_t max_gap_fill = 70) :
        g_(g), path_fixer_(g), skip_unfixed_(skip_unfixed), max_gap_fill_(max_gap_fill)
    {}

    //TODO replace original method, rename sequence() field in path
    Path<EdgeId> FindDetailedReadPath(const MappingPath<EdgeId> &mapping_path) const {
        if (mapping_path.size() == 0) {
            TRACE("Read unmapped");
            return Path<EdgeId>();
        }

        auto fixed_path = path_fixer_.DeleteSameEdges(mapping_path.path());
        if (fixed_path.size() != mapping_path.size()) {
            TRACE("Some edges were deleted");
        }

        fixed_path = path_fixer_.TryFixPath(fixed_path, max_gap_fill_);

        if (!CheckContiguous(g_, fixed_path.sequence())) {
            TRACE("Could not fix the path!")
            if (skip_unfixed_) {
                TRACE("Read unmapped");
                return Path<EdgeId>();
            } else {
                TRACE("Could not fix the path!")
            }
        }
        return fixed_path;
    }

    std::vector<EdgeId> FindReadPath(const MappingPath<EdgeId> &mapping_path) const {
        return FindDetailedReadPath(mapping_path).sequence();
    }
};

template<class Graph, class Index>
class BasicSequenceMapper: public AbstractSequenceMapper<Graph> {
  using AbstractSequenceMapper<Graph>::g_;

  const Index& index_;

  typedef std::vector<MappingRange> RangeMappings;
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  typedef typename Index::KMer Kmer;
  typedef KmerMapper<Graph> KmerSubs;
  const KmerSubs& kmer_mapper_;
  size_t k_;
  bool optimization_on_;

  bool FindKmer(const Kmer &kmer, size_t kmer_pos, std::vector<EdgeId> &passed,
                RangeMappings& range_mappings) const {
    const auto& position = index_.get(kmer);
    if (position.second == Index::NOT_FOUND)
        return false;
    
    if (passed.empty() || passed.back() != position.first ||
        kmer_pos != range_mappings.back().initial_range.end_pos ||
        position.second + 1 < range_mappings.back().mapped_range.end_pos) {
        passed.push_back(position.first);

        range_mappings.emplace_back(Range(kmer_pos, kmer_pos + 1),
                                    Range(position.second, position.second + 1));
    } else {
        range_mappings.back().initial_range.end_pos = kmer_pos + 1;
        range_mappings.back().mapped_range.end_pos = position.second + 1;
    }

    return true;
  }

  bool TryThread(const Kmer& kmer, size_t kmer_pos, std::vector<EdgeId> &passed,
                 RangeMappings& range_mappings) const {
    EdgeId last_edge = passed.back();
    size_t end_pos = range_mappings.back().mapped_range.end_pos;
    if (end_pos < g_.length(last_edge)) {
      const Sequence &seq = g_.EdgeNucls(last_edge);
      if (seq[end_pos + k_ - 1] == kmer[k_ - 1]) {
        range_mappings.back().initial_range.end_pos++;
        range_mappings.back().mapped_range.end_pos++;
        return true;
      }
    } else {
      VertexId v = g_.EdgeEnd(last_edge);

      if (!optimization_on_)
          if (g_.OutgoingEdgeCount(v) > 1)
              return false;

      for (EdgeId edge : g_.OutgoingEdges(v)) {
        const Sequence &seq = g_.EdgeNucls(edge);
        if (seq[k_ - 1] == kmer[k_ - 1]) {
          passed.push_back(edge);
          range_mappings.emplace_back(Range(kmer_pos, kmer_pos + 1),
                                      Range(0, 1));
          return true;
        }
      }
    }
    return false;
  }

  bool ProcessKmer(const Kmer &kmer, size_t kmer_pos, std::vector<EdgeId> &passed_edges,
                   RangeMappings& range_mapping, bool try_thread) const {
    if (try_thread) {
        if (!TryThread(kmer, kmer_pos, passed_edges, range_mapping)) {
            FindKmer(kmer_mapper_.Substitute(kmer), kmer_pos, passed_edges, range_mapping);
            return false;
        }

        return true;
    }

    if (kmer_mapper_.CanSubstitute(kmer)) {
        FindKmer(kmer_mapper_.Substitute(kmer), kmer_pos, passed_edges, range_mapping);
        return false;
    }

    return FindKmer(kmer, kmer_pos, passed_edges, range_mapping);
  }

 public:
  BasicSequenceMapper(const Graph& g,
                      const Index& index,
                      const KmerSubs& kmer_mapper,
                      bool optimization_on = true) :
      AbstractSequenceMapper<Graph>(g), index_(index),
      kmer_mapper_(kmer_mapper), k_(g.k()+1),
      optimization_on_(optimization_on) { }

  MappingPath<EdgeId> MapSequence(const Sequence &sequence,
                                  bool only_simple = false) const {
    std::vector<EdgeId> passed_edges;
    RangeMappings range_mapping;
    
    if (sequence.size() < k_) {
      return MappingPath<EdgeId>();
    }

    Kmer kmer = sequence.start<Kmer>(k_);
    bool try_thread = false;
    try_thread = ProcessKmer(kmer, 0, passed_edges,
                             range_mapping, try_thread);
    for (size_t i = k_; i < sequence.size(); ++i) {
      kmer <<= sequence[i];
      try_thread = ProcessKmer(kmer, i - k_ + 1, passed_edges,
                               range_mapping, try_thread);
      if (only_simple && passed_edges.size() > 1)
        return MappingPath<EdgeId>();
    }

    return MappingPath<EdgeId>(passed_edges, range_mapping);
  }

  DECL_LOGGER("BasicSequenceMapper");
};

std::shared_ptr<BasicSequenceMapper<Graph, EdgeIndex<Graph>>> MapperInstance(const GraphPack &gp);
std::shared_ptr<BasicSequenceMapper<Graph, EdgeIndex<Graph>>> MapperInstance(const GraphPack &gp,
                                                                             const EdgeIndex<Graph> &index);
} // namespace debruijn_graph
