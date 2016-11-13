//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <cstring>
#include <vector>
#include <algorithm>
#include "sequence/sequence.hpp"
#include "sequence/sequence_tools.hpp"

namespace cap {

namespace utils {

// For some internal maps
struct unordered_map_pair_hash {
    size_t operator()(std::pair<unsigned, size_t> p) const {
        return size_t(p.first) | (p.second << 8);
    }
};

bool compare_pairs_reversed(const std::pair<size_t, size_t> p1,
                            const std::pair<size_t, size_t> p2) {
    if (p1.second != p2.second) {
        return p1.second < p2.second;
    }
    return p1.first < p2.first;
}

}

namespace debug {

const unsigned kShiftValue = 24;
const size_t kMaskValue = (1ull << kShiftValue) - 1;

template<class T, class Graph>
std::string Debug(const Graph *g_, const std::vector<T> &p) {
    std::stringstream ss;
    for (const auto &x : p) {
        ss << g_->str(x) << ";";
    }
    return ss.str();
}

template<class T, class Graph>
std::string Debug(const Graph &g_, const std::vector<T> &p) {
    return Debug(&g_, p);
}

std::string PrintComplexPosition(const size_t pos) {
    std::stringstream ss;
    ss << (pos >> kShiftValue) << ":" << (double(pos & kMaskValue) / (1 << kShiftValue));
    return ss.str();
}

std::string PrintComplexRange(const Range &range) {
    return "[" + PrintComplexPosition(range.start_pos) + ", " +
        PrintComplexPosition(range.end_pos) + "]";
}

std::string PrintComplexRange(const std::pair<size_t, size_t> &range) {
    return "[" + PrintComplexPosition(range.first) + ", " +
        PrintComplexPosition(range.second) + "]";
}


}

/*
 *  Since some time positions do not exactly correspond to the positions
 *  in the genomes. Instead, they are shifted left by 16 bits and, probably,
 *  added with some small number in order to exclude empty ranges.
 */
template<class Graph>
class CoordinatesHandler : public ActionHandler<typename Graph::VertexId,
    typename Graph::EdgeId> {
     public:
      typedef typename Graph::EdgeId EdgeId;
      typedef typename Graph::VertexId VertexId;
      typedef ActionHandler<VertexId, EdgeId> base;
      typedef std::vector<EdgeId> Path;
      typedef unsigned uint;
      typedef std::vector<std::pair<uint, size_t> > PosArray;
      typedef std::vector<std::pair<uint, Range> > RangeArray;
      typedef std::vector<std::pair<size_t, size_t> > Thread;

      const static unsigned kShiftValue = 24;
      const static size_t kMaskValue = (1ull << kShiftValue) - 1;
      const static size_t kNotMaskValue = ~kMaskValue;
      const static size_t kHalfMask = (1ull << (kShiftValue - 1));
      const static size_t kLeftEndMask = (1ull << (kShiftValue - 1)) | (1ull << (kShiftValue - 2));
      const static size_t kRightEndMask = (0ull << (kShiftValue - 1)) | (1ull << (kShiftValue - 2));

      CoordinatesHandler()
          : base("CoordinatesHandler"),
          g_(NULL),
          genome_info_(),
          edge_ranges_(),
          stored_threading_history_(),
          //genome_first_edges_(),
          last_deleted_(),
          pending_add_(),
          is_locked_(false) {
          }

      CoordinatesHandler(
          const std::vector<std::pair<uint, std::vector<Thread>>> &stored_threads)
          : CoordinatesHandler() {
              SetStoredThreads(stored_threads);
          }

      virtual ~CoordinatesHandler() {
      }

      void SetGraph(const Graph *g) {
          VERIFY(g != NULL);

          if (g == g_)
              return;

          if (g_ != NULL)
              UnsetGraph();
          g_ = g;
          g_->AddActionHandler(this);
      }

      void UnsetGraph() {
          if (g_ == NULL)
              return;

          g_->RemoveActionHandler(this);
          g_ = NULL;
          edge_ranges_.clear();
      }

      const Graph *GetGraph() const {
          return g_;
      }

      void AddGenomePath(const uint genome_id, const Path &genome_path) {
          TRACE("AddGenomePath Start");
          VERIFY(g_ != NULL);

          if (genome_path.size() == 0) {
              INFO("Trying to add path of length 0");
              return;
          }
          genome_info_[genome_id].id = genome_id;
          genome_info_[genome_id].first_edge = genome_path[0];
          //genome_first_edges_[genome_path[0]].push_back(genome_id);

          size_t cur_start = 0;
          for (const auto &edge : genome_path) {
              if (edge == EdgeId(0)) {
                  DEBUG("ZERO EDGE!");
                  continue;
              }
              const size_t cur_end = (cur_start + (g_->length(edge) << kShiftValue)) | kHalfMask;

              //DEBUG("edge " << g_->str(edge) << ": " << PrintComplexRange(Range(cur_start, cur_end)));

              edge_ranges_[edge].AddGenomeRange(genome_id, Range(cur_start, cur_end));
              cur_start = cur_end;
          }

          genome_info_[genome_id].sequence_length = cur_start;
          TRACE("AddGenomePath End");
      }

      Sequence ReconstructGenome(const uint genome_id) const {
          const std::vector<EdgeId> genome_path =
              AsMappingPath(genome_id).simple_path();

          std::vector<Sequence> path_sequences;
          for (const auto &e : genome_path)
              path_sequences.push_back(g_->EdgeNucls(e));

          return MergeOverlappingSequences(path_sequences, g_->k());
      }

      PosArray FilterPosArray(const PosArray &old_array,
                              const EdgeId edge) const {
          PosArray result;
          auto edge_data_it = edge_ranges_.find(edge);
          if (edge_data_it == edge_ranges_.end())
              return result;

          for (auto entry : old_array) {
              if (edge_data_it->second.HasForwardLink(entry)) {
                  entry.second = edge_data_it->second.GetForwardPos(entry);
                  result.push_back(entry);
              }
          }

          return result;
      }

      void SetStoredThreads(
          const std::vector<std::pair<uint, std::vector<Thread>>> &threads) {
          for (const auto &entry : threads) {
              stored_threading_history_[entry.first] = entry.second;
          }
      }

      std::vector<std::pair<uint, std::vector<Thread>>> GetStoredThreads() const {
          std::vector<std::pair<uint, std::vector<Thread>>> result;
          for (const auto &entry : stored_threading_history_) {
              result.push_back(entry);
          }
          return result;
      }

      /**
       * Automatically adds conjugate strand!!!
       */
      void StoreGenomeThreadManual(const uint genome_id, const Thread &ladder) {
        stored_threading_history_[2 * genome_id].push_back(PreprocessCoordinates(ladder));
        stored_threading_history_[(2 * genome_id) ^ 1].push_back(PreprocessCoordinates(ConjugateThread(ladder)));
      }

      size_t PreprocessCoordinates(const size_t coord) const {
          return (coord << kShiftValue);
      }

      /*
       * `from` is meant to have needed range data
       */
      bool ProjectPath(const Path &from, const Path &to);

      bool ProjectPath(const Path &from, const Path &to,
                       const PosArray &threads_to_delete);

      void UnrollChanges();

      void LockChanges();

      void ReleaseChanges();

      PosArray GetContiguousThreads(const Path &path) const;

      bool CheckCorrectPathProjection(const Path &from, const Path &to) const;

      size_t GetOriginalPos(const uint genome_id, const size_t new_pos) const;

      //size_t GetNewestPos(const uint genome_id, const size_t old_pos) const;

      // TODO getOrigRange?? (Edge)
      virtual void HandleDelete(EdgeId e);

      virtual void HandleMerge(const vector<EdgeId> &old_edges, EdgeId new_edge);

      virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2);

      virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
                               EdgeId new_edge_2);

      EdgeId FindGenomeFirstEdge(const uint genome_id) const;

      // First range is graph range, second is original sequence range
      std::vector<std::pair<uint, std::pair<Range, Range> > > GetRanges(
          const EdgeId e) const {
          std::vector<std::pair<uint, std::pair<Range, Range> > > res;
          if (g_ == NULL)
              return res;

          const auto edge_data_it = edge_ranges_.find(e);
          VERIFY(edge_data_it != edge_ranges_.end());

          for (const auto &p : edge_data_it->second.GetRanges()) {
              const uint genome_id = p.first;
              Range newest = p.second;
              DEBUG("from " << debug::PrintComplexRange(newest));
              const Range original(GetOriginalPos(genome_id, newest.start_pos),
                                   GetOriginalPos(genome_id, newest.end_pos));

              res.push_back(make_pair(genome_id, make_pair(newest, original)));
          }

          return res;
      }

      std::vector<std::pair<uint, Range> > GetRawRanges(const EdgeId e) const {
          if (g_ == NULL)
              return {};

          const auto edge_data_it = edge_ranges_.find(e);
          VERIFY(edge_data_it != edge_ranges_.end());

          return edge_data_it->second.GetRanges();
      }


      static Range GetPrintableRange(const Range &r) {
          return Range(r.start_pos >> kShiftValue, r.end_pos >> kShiftValue);
      }

      PosArray GetEndPosArray(const EdgeId e) const {
          PosArray result;

          const auto edge_data_it = edge_ranges_.find(e);
          if (edge_data_it == edge_ranges_.end())
              return result;

          for (const auto &entry : edge_data_it->second.GetRanges()) {
              result.push_back(make_pair(entry.first, entry.second.end_pos));
          }

          return result;
      }

      bool HasForwardLink(const EdgeId edge, const uint genome_id,
                          const size_t start_pos) const {
          auto edge_it = edge_ranges_.find(edge);
          // VERIFY(edge_it != edge_ranges_.end());
          if (edge_it == edge_ranges_.end())
            return false;

          return edge_it->second.HasForwardLink(make_pair(genome_id, start_pos));
      }

      size_t GetForwardPos(const EdgeId edge, const uint genome_id,
                           const size_t start_pos) const {
          auto edge_it = edge_ranges_.find(edge);
          VERIFY(edge_it != edge_ranges_.end());

          return edge_it->second.GetForwardPos(make_pair(genome_id, start_pos));
      }

      pair<EdgeId, size_t> StepForward(const VertexId v, const uint genome_id,
                                       const size_t pos) const {
          for (EdgeId e : g_->OutgoingEdges(v)) {
              if (HasForwardLink(e, genome_id, pos)) {
                  return make_pair(e, GetForwardPos(e, genome_id, pos));
              }
          }
          return make_pair(EdgeId(0), -1u);
      }

      pair<EdgeId, size_t> StepForwardPos(const EdgeId last_edge, const uint genome_id,
                                          const size_t last_pos) const {
          return StepForwardPos(g_->EdgeEnd(last_edge), genome_id, last_pos);
      }

      size_t GetMultiplicity(const EdgeId edge) const {
          auto edge_it = edge_ranges_.find(edge);
          if (edge_it == edge_ranges_.end())
              return 0;
          return edge_it->second.GetMultiplicity();
      }

      void DebugOutput(const EdgeId e) {
          if (HasEdgeData(e)) {
              DEBUG("edge " << g_->str(e) << " " << edge_ranges_.at(e).DebugOutput());
          } else {
              DEBUG("edge " << g_->str(e) << " empty");
          }
      }

    //todo some usages do not need original pos, optimize if needed
    MappingPath<EdgeId> AsMappingPath(unsigned genome_id) const {
        MappingPath<EdgeId> answer;
        VertexId v = g_->EdgeStart(FindGenomeFirstEdge(genome_id));
        size_t genome_pos = 0;

        while (true) {
            auto step = StepForward(v, genome_id, genome_pos);
            if (step.second == -1u)
                break;
            EdgeId e = step.first;

            size_t next_genome_pos = step.second;

            Range original_pos(
                    GetOriginalPos(genome_id, genome_pos),
                    GetOriginalPos(genome_id, next_genome_pos));

            //todo fix possible troubles with cyclic genomes etc later
            Range original_pos_printable = GetPrintableRange(original_pos);
            Range graph_pos_printable(0, g_->length(e));

            answer.push_back(e, MappingRange(original_pos_printable, graph_pos_printable));

            v = g_->EdgeEnd(e);
            genome_pos = next_genome_pos;
        }
        //todo can we verify total length somehow
        return answer;
    }

    void DumpRanges() {
        std::unordered_map<uint, Path> genome_paths;

        for (const auto &genome_i : genome_info_) {
            const uint genome_id = genome_i.first;

            genome_paths[genome_id] = AsMappingPath(genome_id).
                simple_path();
        }

        StoreGenomeThreads();
        edge_ranges_.clear();
        genome_info_.clear();

        for (const auto &genome_it : genome_paths) {
            const uint genome_id = genome_it.first;

            AddGenomePath(genome_id, genome_it.second);
        }
    }

     private:
      typedef std::unordered_map<std::pair<uint, size_t>, size_t,
              utils::unordered_map_pair_hash> MapT;

      class EdgeData {
       public:
        EdgeData() {
        }

        inline void AddGenomeRange(const uint genome_id, const Range &range) {
            //INFO("add genome range: genome " << int(genome_id) << ": " << range);
            //INFO("initial " << DebugOutput());
            Range extended_range = range;
            auto connected_it = genome_ranges_backward_.find(
                make_pair(genome_id, extended_range.start_pos));
            if (connected_it != genome_ranges_backward_.end()) {
                extended_range.start_pos = connected_it->second;

                genome_ranges_forward_.erase(
                    make_pair(genome_id, extended_range.start_pos));
                genome_ranges_backward_.erase(connected_it);
            }

            connected_it = genome_ranges_forward_.find(
                make_pair(genome_id, extended_range.end_pos));
            if (connected_it != genome_ranges_forward_.end()) {
                extended_range.end_pos = connected_it->second;

                genome_ranges_forward_.erase(connected_it);
                genome_ranges_backward_.erase(
                    make_pair(genome_id, extended_range.end_pos));
            }

            genome_ranges_forward_[make_pair(genome_id, extended_range.start_pos)] =
                extended_range.end_pos;
            genome_ranges_backward_[make_pair(genome_id, extended_range.end_pos)] =
                extended_range.start_pos;

            //INFO("resulting " << DebugOutput());
        }

        inline void operator+=(const EdgeData &other) {
            for (const auto &it : other.genome_ranges_forward_) {
                this->AddGenomeRange(it.first.first, Range(it.first.second, it.second));
            }
        }

        inline std::vector<Range> GetGenomeRanges(const uint genome_id) const {
            std::vector<Range> result;

            for (auto &it : genome_ranges_forward_) {
                if (it.first.first == genome_id) {
                    result.push_back(Range(it.first.second, it.second));
                }
            }

            return result;
        }

        inline std::vector<std::pair<uint, Range> > GetRanges() const {
            std::vector<std::pair<uint, Range> > result;

            for (auto &it : genome_ranges_forward_) {
                result.push_back(make_pair(it.first.first,
                                           Range(it.first.second, it.second)));
            }

            return result;
        }

        inline bool HasForwardLink(const std::pair<uint, size_t> &start_pos) const {
            return genome_ranges_forward_.count(start_pos) > 0;
        }
        inline void DeleteForwardLink(const std::pair<uint, size_t> &pos) {
            const auto it = genome_ranges_forward_.find(pos);
            if (it == genome_ranges_forward_.end())
                return;

            size_t end_pos = it->second;
            genome_ranges_forward_.erase(pos);
            genome_ranges_backward_.erase(make_pair(pos.first, end_pos));
        }
        inline Range PopForwardLink(const std::pair<uint, size_t> &pos) {
            const auto it = genome_ranges_forward_.find(pos);
            if (it == genome_ranges_forward_.end())
                return Range(-1, -1);

            size_t end_pos = it->second;
            genome_ranges_forward_.erase(pos);
            genome_ranges_backward_.erase(make_pair(pos.first, end_pos));

            return Range(pos.second, end_pos);
        }
        inline size_t GetForwardPos(const std::pair<uint, size_t> &start_pos) const {
            auto it = genome_ranges_forward_.find(start_pos);
            VERIFY(it != genome_ranges_forward_.end());
            return it->second;
        }
        inline size_t GetMultiplicity() const {
            return genome_ranges_forward_.size();
        }

        std::string DebugOutput() const {
            using debug::PrintComplexPosition;

            std::stringstream ss;
            ss << "forward (";
            for (const auto &p : genome_ranges_forward_) {
                ss << "{" << int(p.first.first) << "," << PrintComplexPosition(p.first.second) << "}->" << PrintComplexPosition(p.second) << ";";
            }
            ss << ") backward (";
            for (const auto &p : genome_ranges_backward_) {
                ss << "{" << int(p.first.first) << "," << PrintComplexPosition(p.first.second) << "}->" << PrintComplexPosition(p.second) << ";";
            }
            ss << ")";
            return ss.str();
        }

       private:

        MapT genome_ranges_forward_;
        MapT genome_ranges_backward_;
      };

      struct GenomeInfo {
          uint id;
          EdgeId first_edge;
          size_t sequence_length;
      };

      constexpr static long double EPS = 1e-9;

      pair<size_t, size_t> PreprocessCoordinates(const pair<size_t, size_t>& point) const {
          return make_pair(PreprocessCoordinates(point.first), PreprocessCoordinates(point.second));
      }

      Thread PreprocessCoordinates(const Thread& ladder) const {
          Thread answer;
          for (pair<size_t, size_t> point : ladder) {
              answer.push_back(PreprocessCoordinates(point));
          }
          return answer;
      }

      Thread ConjugateThread(const Thread& ladder) const {
          Thread answer;
          size_t first_length = ladder.back().first;
          size_t second_length = ladder.back().second;
          for (auto it = ladder.rbegin(); it != ladder.rend(); ++it) {
              answer.push_back(make_pair(first_length - it->first,
                                         second_length - it->second));
          }
          return answer;
      }

      void StoreGenomeThreads() {
          std::vector<std::pair<std::pair<uint, Range>, EdgeId> > all_ranges;
          /*
           *  A kind of verification
           */
          /*
             for (auto it = g_->SmartEdgeBegin(); !it.IsEnd(); ++it) {
             if (!HasEdgeData(*it)) {
             INFO("NO data for edge " << g_->str(*it));
             }
             for (const auto &range_data : edge_ranges_.at(*it).GetRanges()) {
             all_ranges.push_back(make_pair(make_pair(range_data.first, range_data.second), *it));
             }
             }
             std::sort(all_ranges.begin(), all_ranges.end());
             for (const auto &e : all_ranges) {
             INFO("genome " << int(e.first.first) << ", " << e.first.second << ": " <<
             g_->str(e.second));
             }
             for (size_t i = 1; i < all_ranges.size(); ++i) {
             if (all_ranges[i].first.first != all_ranges[i - 1].first.first) continue;
             if (all_ranges[i].first.second.start_pos != all_ranges[i - 1].first.second.end_pos) {
             INFO("!!! TORN in genome " << int(all_ranges[i].first.first) << " at position " << all_ranges[i].first.second.start_pos);
             }
             }
             */

          TRACE("StoreGenomeThreads Start");
          VERIFY(g_ != NULL);
          for (auto &genome_it : genome_info_) {
              const uint genome_id = genome_it.first;
              stored_threading_history_[genome_id].push_back(Thread());
              Thread &thread = stored_threading_history_[genome_id].back();

              StoreGenomeThread(genome_id, thread);
          }
          TRACE("StoreGenomeThreads End");
      }


      void StoreGenomeThread(const uint genome_id, Thread &thread);

      RangeArray PopAndUpdateRangesToCopy(const EdgeId edge,
                                          PosArray &delete_positions);

      bool CheckContiguousPath(const Path &path) const;

      size_t FlushZeroSequence(const std::vector<EdgeId> &to_edges,
                       const uint genome_id, const Range &from_range,
                       const bool finalize_range);

      size_t GetNonzeroSplit(const Range &from_range, const size_t taken_length,
        const bool finalize_range);

      size_t GetRightEnd(const bool finalize_range, const Range &range) const;

      void AddPendingRangeToEdge(const EdgeId edge, const uint genome_id, const Range &range);

      void DeleteRangeFromEdge(const EdgeId edge, const uint genome_id, const size_t range_from);

      void DeleteRangeFromEdge(const EdgeId edge, const std::pair<uint, size_t> &pos);

      void FlushPendingRanges();

      bool HasEdgeData(const EdgeId edge) {
          return edge_ranges_.find(edge) != edge_ranges_.end();
      }

      template <class T>
    inline void CleanEdgeData(const T edge) {
        edge_ranges_.erase(edge);
    }

    size_t GetPathLength(const Path &p) const {
        size_t res = 0;
        for (const auto &edge : p) {
            res += g_->length(edge);
        }
        return res;
    }

    size_t CalculatePos(const Range &graph_range, const Range &genome_range,
                        const size_t graph_pos) const {
        const long double len_ratio = (long double)GetPrintableRange(genome_range).size() /
            GetPrintableRange(graph_range).size();

        return ((genome_range.start_pos & kNotMaskValue) | kHalfMask) +
            (size_t(((graph_pos - graph_range.start_pos) >> kShiftValue) * len_ratio) << kShiftValue);
    }

    const Graph *g_;

    std::unordered_map<uint, GenomeInfo> genome_info_;
    std::unordered_map<EdgeId, EdgeData> edge_ranges_;
    std::unordered_map<uint, std::vector<Thread> > stored_threading_history_;
    //std::unordered_map<EdgeId, std::vector<uint> > genome_first_edges_;

    // For deletion unroll
    std::vector<std::pair<EdgeId, std::pair<uint, Range> > > last_deleted_;
    // For suspended adding
    std::vector<std::pair<EdgeId, std::pair<uint, Range> > > pending_add_;

    bool is_locked_;

    DECL_LOGGER("CoordinatesHandler");

};

template <class Graph>
bool CoordinatesHandler<Graph>::ProjectPath(const Path &from, const Path &to,
                                            const PosArray &threads_to_delete) {

    using debug::PrintComplexPosition;
    using debug::PrintComplexRange;
    using debug::Debug;

    TRACE("ProjectPath Start");
    //VERIFY(CheckCorrectPathProjection(from, to));
    //DEBUG("Projecting " << Debug(g_, from) << " to " << Debug(g_, to));


    VERIFY(g_ != NULL);
    const Path &p1 = from,
          &p2 = to;
    size_t l1 = GetPathLength(p1),
           l2 = GetPathLength(p2);
    PosArray cur_delete_positions = threads_to_delete;
    const size_t n_ranges = cur_delete_positions.size();

    std::vector<std::vector<EdgeId> > zero_sequences(n_ranges);

    auto it2 = p2.begin();

    VERIFY(l1 != 0 && l2 != 0);
    VERIFY(it2 != p2.end());

    long double lratio = (long double)l2 / l1;

    size_t cur_2_edge_len = g_->length(*it2);
    for (auto &edge1 : p1) {
        size_t cur_1_edge_len = g_->length(edge1);

        RangeArray genome_ranges_to_copy =
            PopAndUpdateRangesToCopy(edge1, cur_delete_positions);

        if (genome_ranges_to_copy.size() != n_ranges) {
            //ClearChanges();
            DEBUG("FALSE; unrolling");
            return false;
        }

        while (cur_1_edge_len > 0 || edge1 == p1.back()) {
            VERIFY(it2 != p2.end());
            const bool second_edge_is_last = (*it2 == p2.back());
            size_t taken_len_1 = min(cur_1_edge_len,
                                     size_t(ceil(double(cur_2_edge_len / lratio - EPS))));
            if (second_edge_is_last) {
                taken_len_1 = cur_1_edge_len;
            }
            const size_t taken_len_2 = min(cur_2_edge_len,
                                           size_t(ceil(double(taken_len_1 * lratio - EPS))));
            const long double edge_1_percentage = (cur_1_edge_len == 0) ? 0 :
                (long double)taken_len_1 / cur_1_edge_len;

            const bool finalize_range = (edge1 == p1.back() && second_edge_is_last)
                        || (edge1 != p1.back() && taken_len_1 == cur_1_edge_len);

            for (size_t i = 0; i < n_ranges; ++i) {
                auto &ranges = genome_ranges_to_copy[i];
                Range &range = ranges.second;
                const size_t range_size = GetPrintableRange(range).size();
                const size_t taken_length =
                    size_t(ceil(double(range_size * edge_1_percentage - EPS)));

                if (taken_length == 0) {
                    zero_sequences[i].push_back(*it2);

                    DEBUG("taken length is zero");
                    continue;
                }

                if (!zero_sequences[i].empty()) {
                    range.start_pos = FlushZeroSequence(zero_sequences[i],
                            ranges.first, range, false);//, taken_length == 0 && finalize_range);
                    zero_sequences[i].clear();
                }

                const size_t split_pos = GetNonzeroSplit(
                    range, taken_length, finalize_range);
                const Range range_to_add(range.start_pos, split_pos);

                /*
                DEBUG("DEBUG: " << PrintComplexPosition(range.start_pos) << " " << taken_length << " "
                      << PrintComplexPosition(split_pos));
                DEBUG("  Proj " << g_->str(edge1) << " -> " << g_->str(*it2) << ": "
                      << int(ranges.first) << ": " << PrintComplexRange(ranges.second) << "->" << PrintComplexRange(range_to_add));
                      */

                AddPendingRangeToEdge(*it2, ranges.first, range_to_add);
                range.start_pos = split_pos;
            }

            cur_1_edge_len -= taken_len_1;
            cur_2_edge_len -= taken_len_2;

            if (edge1 == p1.back()) {
                ++it2;
                if (it2 == p2.end()) {
                    break;
                }
                cur_2_edge_len = g_->length(*it2);
            } else if (cur_2_edge_len == 0) {
                if (*it2 != p2.back()) {
                    ++it2;
                    cur_2_edge_len = g_->length(*it2);
                }
            }
        }


        // Check that all ranges were moved completely
        for (size_t i = 0; i < n_ranges; ++i) {
            const uint genome_id = genome_ranges_to_copy[i].first;
            Range &range = genome_ranges_to_copy[i].second;

            if (!zero_sequences[i].empty()) {
                range.start_pos = FlushZeroSequence(zero_sequences[i],
                        genome_id, range, true);//, taken_length == 0 && finalize_range);
                zero_sequences[i].clear();
            }

            VERIFY(range.start_pos == range.end_pos);
        }
    }

    FlushPendingRanges();

    TRACE("ProjectPath End");

    return true;
}

template <class Graph>
bool CoordinatesHandler<Graph>::ProjectPath(
    const Path &from, const Path &to) {
    PosArray all_positions =
        GetContiguousThreads(from);
    return ProjectPath(from, to, all_positions);
}

template <class Graph>
typename CoordinatesHandler<Graph>::RangeArray
CoordinatesHandler<Graph>::PopAndUpdateRangesToCopy(
    const EdgeId edge,
    PosArray &delete_positions) {
    auto edge_data_it = edge_ranges_.find(edge);
    if (edge_data_it == edge_ranges_.end()) {
        INFO("trying to get " << delete_positions.size() << " positions from empty!!!");
        return {};
    }
    //VERIFY(edge_data_it != edge_ranges_.end());
    auto &edge_data = edge_data_it->second;

    RangeArray genome_ranges_to_copy;
    for (auto &del_pos : delete_positions) {
        if (edge_data.HasForwardLink(del_pos)) {
            const Range range_to_copy(del_pos.second,
                                      edge_data.GetForwardPos(del_pos));
            DeleteRangeFromEdge(edge, del_pos);

            genome_ranges_to_copy.push_back(
                make_pair(del_pos.first, range_to_copy));
            del_pos.second = range_to_copy.end_pos;
        }
    }
    //if (edge_data.GetMultiplicity() == 0)
    //  CleanEdgeData(edge_data_it);

    return genome_ranges_to_copy;
}

template <class Graph>
bool CoordinatesHandler<Graph>::CheckCorrectPathProjection(
    const Path &from, const Path &to) const {
    if (from.size() == 0 || to.size() == 0) return false;
    if (g_->EdgeStart(from[0]) != g_->EdgeStart(to[0])) return false;
    if (g_->EdgeEnd(from.back()) != g_->EdgeEnd(to.back())) return false;
    if (!CheckContiguousPath(from) || !CheckContiguousPath(to)) return false;
    return true;
}

template<class Graph>
size_t CoordinatesHandler<Graph>::FlushZeroSequence(const std::vector<EdgeId> &to_edges,
                       const uint genome_id, const Range &from_range,
                       const bool finalize_range) {
    TRACE("FlushZeroSequence " << debug::Debug(g_, to_edges) << " " << genome_id << " " << debug::PrintComplexRange(from_range) << " " << finalize_range);
    const size_t N = to_edges.size();
    // the least power of 2 that is not smaller than N
    uint K = 0;
    while ((1u << K) < N)
        ++K;

    const size_t right_end = GetRightEnd(finalize_range, from_range);
    VERIFY((right_end & kNotMaskValue) == (from_range.start_pos & kNotMaskValue));

    const size_t dx = right_end - from_range.start_pos;
    const size_t step = dx >> K;
    VERIFY(step != 0);

    size_t cur_pos = from_range.start_pos;

    for (size_t i = 0; i < N; ++i) {
        size_t next_pos = cur_pos + step;
        if (i + 1 == N)
            next_pos = right_end;

        AddPendingRangeToEdge(to_edges[i], genome_id, Range(cur_pos, next_pos));
        cur_pos = next_pos;
    }

    return right_end;
}

template<class Graph>
size_t CoordinatesHandler<Graph>::GetNonzeroSplit(const Range &from_range, const size_t taken_length,
        const bool finalize_range) {
    size_t split_pos = 0;

    if (finalize_range)
        split_pos = from_range.end_pos;
    else {
        // just base
        split_pos = (from_range.start_pos & kNotMaskValue) + (taken_length << kShiftValue);
        if ((from_range.end_pos & kNotMaskValue) != split_pos) {
            split_pos |= kHalfMask;
        } else {
            split_pos |= (from_range.end_pos & kMaskValue) >> 1;
        }
    }

    return split_pos;
}

template<class Graph>
size_t CoordinatesHandler<Graph>::GetRightEnd(const bool finalize_range, const Range &range) const {
    size_t right_end = 0;

    if (finalize_range)
        right_end = range.end_pos;
    else if ((range.start_pos & kMaskValue) == kRightEndMask)
        right_end = (range.start_pos & kNotMaskValue) | kLeftEndMask;
    else {
        const size_t was = range.start_pos & kMaskValue;
        const size_t one_minus_was = (1ull << kShiftValue) - was;
        const size_t new_right = (1ull << kShiftValue) - (one_minus_was >> 1);
        right_end = (range.start_pos & kNotMaskValue) | new_right;
    }

    return right_end;
}

template<class Graph>
void CoordinatesHandler<Graph>::AddPendingRangeToEdge(const EdgeId edge,
        const uint genome_id, const Range &range) {
    pending_add_.push_back(make_pair(edge, make_pair(genome_id, range)));
}

template<class Graph>
void CoordinatesHandler<Graph>::DeleteRangeFromEdge(const EdgeId edge,
        const uint genome_id, const size_t range_from) {
    DeleteRangeFromEdge(edge, make_pair(genome_id, range_from));
}
template<class Graph>
void CoordinatesHandler<Graph>::DeleteRangeFromEdge(const EdgeId edge,
        const std::pair<uint, size_t> &pos) {
    const Range deleted_range =
        edge_ranges_.at(edge).PopForwardLink(pos);
    if (edge_ranges_.at(edge).GetMultiplicity() == 0) {
        CleanEdgeData(edge);
    }

    last_deleted_.push_back(make_pair(edge,
                make_pair(pos.first, deleted_range)));
}

template<class Graph>
void CoordinatesHandler<Graph>::FlushPendingRanges() {
    if (is_locked_)
        return;

    for (const auto &e : pending_add_) {
        edge_ranges_[e.first].AddGenomeRange(e.second.first, e.second.second);
    }
    last_deleted_.clear();
    pending_add_.clear();
}

template<class Graph>
void CoordinatesHandler<Graph>::UnrollChanges() {
    for (const auto &e : last_deleted_) {
        edge_ranges_[e.first].AddGenomeRange(e.second.first, e.second.second);
    }
    pending_add_.clear();
    last_deleted_.clear();
}

template<class Graph>
void CoordinatesHandler<Graph>::ReleaseChanges() {
    is_locked_ = false;
    FlushPendingRanges();
}
template<class Graph>
void CoordinatesHandler<Graph>::LockChanges() {
    is_locked_ = true;
}

template <class Graph>
bool CoordinatesHandler<Graph>::CheckContiguousPath(
    const Path &path) const {
    for (size_t i = 1; i < path.size(); ++i) {
        if (g_->EdgeEnd(path[i - 1]) != g_->EdgeStart(path[i]))
            return false;
    }
    return true;
}

template <class Graph>
typename CoordinatesHandler<Graph>::PosArray
CoordinatesHandler<Graph>::GetContiguousThreads(const Path &path) const {
    PosArray result;
    PosArray cur_pos;

    const auto first_edge_it = edge_ranges_.find(path[0]);
    if (first_edge_it == edge_ranges_.end())
        return result;

    for (const auto &entry : first_edge_it->second.GetRanges()) {
        result.push_back(make_pair(entry.first, entry.second.start_pos));
        cur_pos.push_back(make_pair(entry.first, entry.second.end_pos));
    }

    for (size_t path_i = 1; path_i < path.size(); ++path_i) {
        const auto edge_data_it = edge_ranges_.find(path[path_i]);
        if (edge_data_it == edge_ranges_.end())
            return {};

        for (ssize_t i = 0; i < (ssize_t)cur_pos.size(); ++i) {
            if (edge_data_it->second.HasForwardLink(cur_pos[i])) {
                cur_pos[i].second = edge_data_it->second.GetForwardPos(cur_pos[i]);
            } else {
                std::swap(cur_pos[i], cur_pos.back());
                std::swap(result[i], result.back());
                cur_pos.pop_back();
                result.pop_back();
                i--;
            }
        }
    }

    return result;
}

template <class Graph>
void CoordinatesHandler<Graph>::StoreGenomeThread(
    const uint genome_id, Thread &thread) {

    TRACE("StoreGenomeThread Start");

    size_t graph_pos = 0,
           genome_pos = 0,
           genome_length = genome_info_[genome_id].sequence_length;

    std::pair<uint, size_t> cur_pos = make_pair(genome_id, genome_pos);
    EdgeId cur_edge = genome_info_[genome_id].first_edge;
    //INFO("searching for first edge " << cur_edge << "of genome " << int(genome_id));
    if (!HasEdgeData(cur_edge)) {
        cur_edge = FindGenomeFirstEdge(genome_id);
    }

    do {
        if (cur_edge == EdgeId(0)) {
            INFO("Could not thread genome path! genome_id=" << int(genome_id));
            return;
        }

        thread.push_back(make_pair(graph_pos, genome_pos));
        graph_pos += g_->length(cur_edge) << kShiftValue;
        VERIFY(HasEdgeData(cur_edge));
        genome_pos = edge_ranges_.at(cur_edge).GetForwardPos(cur_pos);
        cur_pos.second = genome_pos;

        const VertexId v = g_->EdgeEnd(cur_edge);

        //DEBUG("current edge " << g_->str(cur_edge) << ", outgoing count " << g_->OutgoingEdgeCount(v));
        cur_edge = EdgeId(0);
        for (const auto &out_edge : g_->OutgoingEdges(v)) {
            //DEBUG("considering edge " << g_->str(out_edge) << " at position (seq) " << genome_pos);

            auto edge_info_it = edge_ranges_.find(out_edge);
            if (edge_info_it == edge_ranges_.end())
                continue;

            //TRACE("!");

            if (edge_info_it->second.HasForwardLink(cur_pos)) {
                cur_edge = out_edge;
                break;
            }
        }
    } while (genome_pos != genome_length);

    thread.push_back(make_pair(graph_pos, genome_pos));
    TRACE("StoreGenomeThread End");
}

template <class Graph>
typename CoordinatesHandler<Graph>::EdgeId
CoordinatesHandler<Graph>::FindGenomeFirstEdge(const uint genome_id) const {
    std::pair<uint, size_t> pos_in_question(genome_id, 0);

    for (auto it = g_->SmartEdgeBegin(); !it.IsEnd(); ++it) {
        auto range_it = edge_ranges_.find(*it);
        if (range_it == edge_ranges_.end())
            continue;
        if (range_it->second.HasForwardLink(pos_in_question)) {
            return *it;
        }
    }

    // remember first edge and update it
    VERIFY_MSG(false, "Could not find start of the sequence in graph");
    return EdgeId(0);
}

template <class Graph>
size_t CoordinatesHandler<Graph>::GetOriginalPos(
    const uint genome_id, const size_t new_pos) const {

    // No refinement has been done
    if (stored_threading_history_.size() == 0)
        return new_pos;

    size_t cur_pos = new_pos;

    const auto history_it = stored_threading_history_.find(genome_id);
    VERIFY(history_it != stored_threading_history_.end());
    const std::vector<Thread> &history = history_it->second;

    for (auto thread_it = history.rbegin(), E = history.rend();
         thread_it != E; ++thread_it) {
        // Verify thread sort order?

        VERIFY(thread_it->size() > 0);

        // Kmers can have different lengths so going from larger kmers to smaller
        // implies shorting of thread length what may lead to "range-overflow"
        if (cur_pos > thread_it->back().first)
            cur_pos = thread_it->back().first;

        auto found_it = std::lower_bound(thread_it->begin(), thread_it->end(),
                                         make_pair(cur_pos, size_t(0)));

        DEBUG("Searching for pos " << debug::PrintComplexPosition(cur_pos)
                << "in thread of " << debug::PrintComplexRange(thread_it->front())
                << " - " << debug::PrintComplexRange(thread_it->back()));
        VERIFY(found_it != thread_it->end());
        if (cur_pos == found_it->first) {
            cur_pos = found_it->second;
            continue;
        }
        VERIFY(found_it != thread_it->begin());

        Range graph_range(0, 0);
        Range genome_range(0, 0);

        graph_range.end_pos = found_it->first;
        genome_range.end_pos = found_it->second;
        --found_it;
        graph_range.start_pos = found_it->first;
        genome_range.start_pos = found_it->second;

        DEBUG("from ranges " << debug::PrintComplexRange(graph_range) <<
                " and " << debug::PrintComplexRange(genome_range) <<
                " in search of " << debug::PrintComplexPosition(cur_pos));
        cur_pos = CalculatePos(graph_range, genome_range, cur_pos);
        DEBUG("gettin' " << debug::PrintComplexPosition(cur_pos));
    }

    return cur_pos;
}

/*
template<class Graph>
size_t CoordinatesHandler<Graph>::GetNewestPos(
    const uint genome_id, const size_t old_pos) const {
    if (stored_threading_history_.size() == 0)
        return old_pos;

    const auto history_it = stored_threading_history_.find(genome_id);
    VERIFY(history_it != stored_threading_history_.end());
    const std::vector<Thread> &history = history_it->second;
    const Thread &latest = history.back();

    // Kmers can have different lengths so going from larger kmers to smaller
    // implies shorting of thread length what may lead to "range-overflow"
    size_t search_pos = old_pos;
    if (search_pos > latest.back().second)
        search_pos = latest.back().second;

    auto found_it = std::lower_bound(latest.begin(), latest.end(),
                                     make_pair(size_t(0), search_pos), utils::compare_pairs_reversed);

    VERIFY(found_it != latest.end());
    if (search_pos == found_it->second)
        return found_it->first;

    VERIFY(found_it != latest.begin());

    Range graph_range(0, 0);
    Range genome_range(0, 0);

    graph_range.end_pos = found_it->second;
    genome_range.end_pos = found_it->first;
    --found_it;
    graph_range.start_pos = found_it->second;
    genome_range.start_pos = found_it->first;

    DEBUG("from ranges " << debug::PrintComplexRange(graph_range) <<
            " and " << debug::PrintComplexRange(genome_range) <<
            " in search of " << debug::PrintComplexPosition(search_pos));
    const size_t result_pos = CalculatePos(graph_range, genome_range, search_pos);
    DEBUG("gettin' " << debug::PrintComplexPosition(result_pos));
    return result_pos;
}
*/

/*
 *  Some handling methods (description in omni_utils.hpp)
 */
template <class Graph>
void CoordinatesHandler<Graph>::HandleDelete(EdgeId e) {
    if (HasEdgeData(e)) {
        INFO("edge " << g_->str(e) << " " << edge_ranges_[e].DebugOutput());
    }
    VERIFY(!HasEdgeData(e));
    CleanEdgeData(e);
}

template <class Graph>
void CoordinatesHandler<Graph>::HandleMerge(const vector<EdgeId> &old_edges, EdgeId new_edge) {
    //DEBUG("HandleMerge : " << debug::Debug(g_, old_edges) << " -> " << g_->str(new_edge));
    for (const auto &edge : old_edges) {
        if (HasEdgeData(edge)) {
            //DEBUG("edge " << g_->str(edge) << " " << edge_ranges_[edge].DebugOutput());
            edge_ranges_[new_edge] += edge_ranges_[edge];
            CleanEdgeData(edge);
        }
    }
}

template <class Graph>
void CoordinatesHandler<Graph>::HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
    //DEBUG("HandleGlue : " << g_->str(new_edge) << " <- " << g_->str(edge1) << " + " << g_->str(edge2));
    if (HasEdgeData(edge1)) {
        //DEBUG("edge " << g_->str(edge1) << " " << edge_ranges_[edge1].DebugOutput());
        edge_ranges_[new_edge] += edge_ranges_[edge1];
        CleanEdgeData(edge1);
    }
    if (HasEdgeData(edge2)) {
        //DEBUG("edge " << g_->str(edge2) << " " << edge_ranges_[edge2].DebugOutput());
        edge_ranges_[new_edge] += edge_ranges_[edge2];
        CleanEdgeData(edge2);
    }
}

template <class Graph>
void CoordinatesHandler<Graph>::HandleSplit(EdgeId old_edge, EdgeId /* new_edge_1 */,
                                            EdgeId /* new_edge_2 */) {
    //DEBUG("HandleSplit " << g_->str(old_edge) << " -> " << g_->str(new_edge_1) << " + " << g_->str(new_edge_2));
    VERIFY(!HasEdgeData(old_edge));
    /*
       const std::vector<std::pair<uint, Range> > old_ranges =
       edge_ranges[old_edge].GetRanges();
       for (const auto &range : old_ranges) {

       }
       */
}


}
