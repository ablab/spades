#pragma once

#include <cstring>
#include <vector>
#include <algorithm>
#include "omni/omni_utils.hpp"
#include "sequence/sequence.hpp"
#include "sequence/sequence_tools.hpp"

namespace cap {

namespace utils {

// For some internal maps
struct unordered_map_pair_hash {
  size_t operator()(std::pair<unsigned char, size_t> p) const {
    //return std::hash<unsigned char>()(p.first) |
        //(std::hash<size_t>()(p.second) << 8);
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

template<class Graph>
class CoordinatesHandler : public ActionHandler<typename Graph::VertexId,
                                                  typename Graph::EdgeId>,
                           private boost::noncopyable {
 public:
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  typedef ActionHandler<VertexId, EdgeId> base;
  typedef std::vector<EdgeId> Path;
  typedef unsigned char uchar;
  typedef std::vector<std::pair<uchar, size_t> > PosArray;

  CoordinatesHandler()
      : base("CoordinatesHandler"),
        g_(NULL),
        genome_info_(),
        edge_ranges_(),
        stored_threading_history_(),
        genome_first_edges_(),
        cur_genome_threads_() {
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

  void AddGenomePath(const uchar genome_id, const Path &genome_path) {
    TRACE("AddGenomePath Start");
    VERIFY(g_ != NULL);

    if (genome_path.size() == 0) {
      INFO("Trying to add path of length 0");
      return;
    }
    genome_info_[genome_id].id = genome_id;
    genome_info_[genome_id].first_edge = genome_path[0];
    genome_first_edges_[genome_path[0]].push_back(genome_id);

    size_t cur_start = 0;
    for (const auto &edge : genome_path) {
      if (edge == EdgeId(0)) {
        DEBUG("ZERO EDGE!");
        continue;
      }
      size_t cur_end = cur_start + g_->length(edge);

      DEBUG("edge " << g_->str(edge) << ": " << Range(cur_start, cur_end));

      edge_ranges_[edge].AddGenomeRange(genome_id, Range(cur_start, cur_end));
      cur_start = cur_end;
    }

    genome_info_[genome_id].sequence_length = cur_start;
    SetGenomeThread(genome_id, genome_path);
    TRACE("AddGenomePath End");
  }

  Sequence ReconstructGenome(const uchar genome_id) const {
    VERIFY(cur_genome_threads_.size() > genome_id);

    std::vector<Sequence> path_sequences;
    path_sequences.reserve(cur_genome_threads_[genome_id]->size());
    for (const auto &e : *cur_genome_threads_[genome_id])
      path_sequences.push_back(g_->EdgeNucls(e));

    return MergeOverlappingSequences(path_sequences, g_->k());
  }

  const Path &GetGenomePath(const uchar genome_id) const {
    VERIFY(cur_genome_threads_.size() > genome_id);

    return *cur_genome_threads_[genome_id];
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

  void StoreGenomeThreads() {
    std::vector<std::pair<std::pair<uchar, Range>, EdgeId> > all_ranges;
    /*
     *  A kind of verification
     */
    for (auto it = g_->SmartEdgeBegin(); !it.IsEnd(); ++it) {
      for (const auto &range_data : edge_ranges_[*it].GetRanges()) {
        all_ranges.push_back(make_pair(make_pair(range_data.first, range_data.second), *it));
      }
    }
    std::sort(all_ranges.begin(), all_ranges.end());
    /*
    for (const auto &e : all_ranges) {
      INFO("genome " << int(e.first.first) << ", " << e.first.second << ": " <<
          g_->str(e.second));
    }
    */
    for (size_t i = 1; i < all_ranges.size(); ++i) {
      if (all_ranges[i].first.first != all_ranges[i - 1].first.first) continue;
      if (all_ranges[i].first.second.start_pos != all_ranges[i - 1].first.second.end_pos) {
        INFO("!!! TORN in genome " << int(all_ranges[i].first.first) << " at position " << all_ranges[i].first.second.start_pos);
      }
    }

    TRACE("StoreGenomeThreads Start");
    VERIFY(g_ != NULL);
    for (auto &genome_it : genome_info_) {
      const uchar genome_id = genome_it.first;
      stored_threading_history_[genome_id].push_back(Thread());
      Thread &thread = stored_threading_history_[genome_id].back();

      StoreGenomeThread(genome_id, thread);
    }
    TRACE("StoreGenomeThreads End");
  }

  /*
   * `from` is meant to have needed range data
   */
  void ProjectPath(const Path &from, const Path &to); 
  void ProjectPath(const Path &from, const Path &to,
      const std::vector<std::pair<uchar, size_t> > &threads_to_delete);

  std::vector<std::pair<uchar, size_t> > GetContiguousThreads(const Path &path) const;
  bool CheckCorrectPathProjection(const Path &from, const Path &to) const;
  size_t GetOriginalPos(const uchar genome_id, const size_t new_pos) const;
  size_t GetNewestPos(const uchar genome_id, const size_t old_pos) const;
  // TODO getOrigRange?? (Edge)
  virtual void HandleDelete(EdgeId e);
  virtual void HandleMerge(const vector<EdgeId> &old_edges, EdgeId new_edge);
  virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2);
  virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
                                            EdgeId new_edge_2);

  // First range is graph range, second is original sequence range
  std::vector<std::pair<uchar, std::pair<Range, Range> > > GetRanges(
      const EdgeId e) const {
    std::vector<std::pair<uchar, std::pair<Range, Range> > > res;
    if (g_ == NULL)
      return res;

    const auto edge_data_it = edge_ranges_.find(e);
    VERIFY(edge_data_it != edge_ranges_.end());

    for (const auto &p : edge_data_it->second.GetRanges()) {
      const uchar genome_id = p.first;
      DEBUG("from " << p.second);
      Range newest(GetNewestPos(genome_id, p.second.start_pos),
                     GetNewestPos(genome_id, p.second.end_pos));
      DEBUG("from " << newest);
      const Range original(GetOriginalPos(genome_id, newest.start_pos),
                             GetOriginalPos(genome_id, newest.end_pos));
      res.push_back(make_pair(genome_id, make_pair(newest, original)));
    }

    return res;
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

  bool HasForwardLink(const EdgeId edge, const uchar genome_id,
                      const size_t start_pos) const {
    auto edge_it = edge_ranges_.find(edge);
    VERIFY(edge_it != edge_ranges_.end());

    return edge_it->second.HasForwardLink(make_pair(genome_id, start_pos));
  }

  size_t GetForwardPos(const EdgeId edge, const uchar genome_id,
                      const size_t start_pos) const {
    auto edge_it = edge_ranges_.find(edge);
    VERIFY(edge_it != edge_ranges_.end());

    return edge_it->second.GetForwardPos(make_pair(genome_id, start_pos));
  }

  size_t GetMultiplicity(const EdgeId edge) const {
    auto edge_it = edge_ranges_.find(edge);
    if (edge_it == edge_ranges_.end())
      return 0;
    return edge_it->second.GetMultiplicity();
  }

  void DebugOutput(const EdgeId e) {
    INFO("edge " << g_->str(e) << " " << edge_ranges_[e].DebugOutput());
  }

 private:
  typedef std::unordered_map<std::pair<uchar, size_t>, size_t,
                             utils::unordered_map_pair_hash> MapT;
  typedef std::vector<std::pair<size_t, size_t> > Thread;

  class EdgeData {
   public:
    EdgeData() {
    }

    inline void AddGenomeRange(const uchar genome_id, const Range &range) {
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

    inline std::vector<Range> GetGenomeRanges(const uchar genome_id) const {
      std::vector<Range> result;

      for (auto &it : genome_ranges_forward_) {
        if (it.first.first == genome_id) {
          result.push_back(Range(it.first.second, it.second));
        }
      }

      return result;
    }

    inline std::vector<std::pair<uchar, Range> > GetRanges() const {
      std::vector<std::pair<uchar, Range> > result;
      
      for (auto &it : genome_ranges_forward_) {
        result.push_back(make_pair(it.first.first,
                                   Range(it.first.second, it.second)));
      }

      return result;
    }

    inline bool HasForwardLink(const std::pair<uchar, size_t> &start_pos) const {
      return genome_ranges_forward_.count(start_pos) > 0;
    }
    inline void DeleteForwardLink(const std::pair<uchar, size_t> &pos) {
      const auto it = genome_ranges_forward_.find(pos);
      if (it == genome_ranges_forward_.end())
        return;

      size_t end_pos = it->second;
      genome_ranges_forward_.erase(pos);
      genome_ranges_backward_.erase(make_pair(pos.first, end_pos));
    }
    inline size_t GetForwardPos(const std::pair<uchar, size_t> &start_pos) const {
      auto it = genome_ranges_forward_.find(start_pos);
      VERIFY(it != genome_ranges_forward_.end());
      return it->second;
    }
    inline size_t GetMultiplicity() const {
      return genome_ranges_forward_.size();
    }

    std::string DebugOutput() const {
      std::stringstream ss;
      ss << "forward (";
      for (const auto &p : genome_ranges_forward_) {
        ss << "{" << int(p.first.first) << "," << p.first.second << "}->" << p.second << ";";
      }
      ss << ") backward (";
      for (const auto &p : genome_ranges_backward_) {
        ss << "{" << int(p.first.first) << "," << p.first.second << "}->" << p.second << ";";
      }
      ss << ")";
      return ss.str();
    }

   private:

    MapT genome_ranges_forward_;
    MapT genome_ranges_backward_;
  };

  struct GenomeInfo {
    uchar id;
    EdgeId first_edge;
    size_t sequence_length;
  };

  constexpr static long double EPS = 1e-9;

  void StoreGenomeThread(const uchar genome_id, Thread &thread);
  EdgeId FindGenomeFirstEdge(const uchar genome_id) const;
  std::vector<std::pair<uchar, Range> > PopAndUpdateRangesToCopy(
      const EdgeId edge,
      std::vector<std::pair<uchar, size_t> > &delete_positions);
  bool CheckContiguousPath(const Path &path) const;

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
    long double len_ratio = (long double)genome_range.size() /
                                         graph_range.size();

    return size_t(genome_range.start_pos +
                  (graph_pos - graph_range.start_pos) * len_ratio);
  }

  void SetGenomeThread(const uchar genome_id, const Path &thread) {
    if (cur_genome_threads_.size() <= genome_id)
      cur_genome_threads_.resize(genome_id + 1);
    cur_genome_threads_[genome_id] = std::make_shared<Path>(thread);
  }

  

  const Graph *g_;
  
  std::unordered_map<uchar, GenomeInfo> genome_info_;
  std::unordered_map<EdgeId, EdgeData> edge_ranges_;
  std::unordered_map<uchar, std::vector<Thread> > stored_threading_history_;
  std::unordered_map<EdgeId, std::vector<uchar> > genome_first_edges_;
  std::vector<std::shared_ptr<Path> > cur_genome_threads_;

  DECL_LOGGER("CoordinatesHandler")
    ;
  template<class T>
  std::string Debug(const std::vector<T> &p) {
    std::stringstream ss;
    for (const auto &x : p) {
      ss << g_->str(x) << ";";
    }
    return ss.str();
  }
};

template <class Graph>
void CoordinatesHandler<Graph>::ProjectPath(const Path &from, const Path &to,
    const std::vector<std::pair<uchar, size_t> > &threads_to_delete) {

  TRACE("ProjectPath Start");
  //VERIFY(CheckCorrectPathProjection(from, to));
  //INFO("Projecting " << Debug(from) << " to " << Debug(to));


  VERIFY(g_ != NULL);
  const Path &p1 = from,
             &p2 = to;
  size_t l1 = GetPathLength(p1),
         l2 = GetPathLength(p2);
  std::vector<std::pair<uchar, size_t> > cur_delete_positions =
      threads_to_delete;
  // For suspended adding
  std::vector<std::pair<EdgeId, std::pair<uchar, Range> > > adding_ranges;

  auto it2 = p2.begin();

  VERIFY(l1 != 0 && l2 != 0);
  VERIFY(it2 != p2.end());

  long double lratio = (long double)l2 / l1;

  size_t cur_2_edge_len = g_->length(*it2);
  for (auto &edge1 : p1) {
    size_t cur_1_edge_len = g_->length(edge1);

    std::vector<std::pair<uchar, Range> > genome_ranges_to_copy =
      PopAndUpdateRangesToCopy(edge1, cur_delete_positions);
    //VERIFY(genome_ranges_to_copy.size() == threads_to_delete.size());

    while (cur_1_edge_len > 0 || edge1 == p1.back()) {
      VERIFY(it2 != p2.end());
      size_t taken_len_1 = min(cur_1_edge_len,
                             size_t(ceil(cur_2_edge_len / lratio - EPS)));
      if (*it2 == p2.back()) {
        taken_len_1 = cur_1_edge_len;
      }
      const size_t taken_len_2 = min(cur_2_edge_len,
                             size_t(ceil(taken_len_1 * lratio - EPS)));
      const long double edge_1_percentage = (cur_1_edge_len == 0) ? 0 :
          (long double)taken_len_1 / cur_1_edge_len;
      
      for (auto &ranges : genome_ranges_to_copy) {
        const size_t taken_length =
            size_t(ceil(ranges.second.size() * edge_1_percentage - EPS));
        //VERIFY(taken_length > 0);
        const Range range_to_add(ranges.second.start_pos,
            ranges.second.start_pos + taken_length);

        //INFO("  Proj " << g_->str(edge1) << " -> " << g_->str(*it2) << ": "
        //    << int(ranges.first) << ": " << ranges.second << "->" << range_to_add);
        //INFO("DEBUG: " << range_to_add << " from range " << ranges.second);
        //edge_ranges_[*it2].AddGenomeRange(ranges.first, range_to_add);
        adding_ranges.push_back(make_pair(*it2, make_pair(ranges.first, range_to_add)));
        ranges.second.start_pos += taken_length;
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
  }

  for (const auto &e : adding_ranges) {
    edge_ranges_[e.first].AddGenomeRange(e.second.first, e.second.second);
  }

  TRACE("ProjectPath End");
}
template <class Graph>
void CoordinatesHandler<Graph>::ProjectPath(
    const Path &from, const Path &to) {
  std::vector<std::pair<uchar, size_t> > all_positions =
      GetContiguousThreads(from);
  ProjectPath(from, to, all_positions);
}

template <class Graph>
std::vector<std::pair<unsigned char, Range> >
CoordinatesHandler<Graph>::PopAndUpdateRangesToCopy(
    const EdgeId edge,
    std::vector<std::pair<uchar, size_t> > &delete_positions) {
  auto edge_data_it = edge_ranges_.find(edge);
  if (edge_data_it == edge_ranges_.end()) {
    INFO("trying to get " << delete_positions.size() << " positions from empty!!!");
  }
  VERIFY(edge_data_it != edge_ranges_.end());
  auto &edge_data = edge_data_it->second;

  std::vector<std::pair<uchar, Range> > genome_ranges_to_copy;
  for (auto &del_pos : delete_positions) {
    if (edge_data.HasForwardLink(del_pos)) {
      const Range range_to_copy(del_pos.second,
          edge_data.GetForwardPos(del_pos));
      edge_data.DeleteForwardLink(del_pos);

      genome_ranges_to_copy.push_back(
          make_pair(del_pos.first, range_to_copy));
      del_pos.second = range_to_copy.end_pos;
    }
  }
  if (edge_data.GetMultiplicity() == 0)
    CleanEdgeData(edge_data_it);

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

    for (size_t i = 0; i < cur_pos.size(); ++i) {
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
    const uchar genome_id, Thread &thread) {

  TRACE("StoreGenomeThread Start");

  size_t graph_pos = 0,
         genome_pos = 0,
         genome_length = genome_info_[genome_id].sequence_length;

  cur_genome_threads_[genome_id] = std::make_shared<Path>();
  std::pair<uchar, size_t> cur_pos = make_pair(genome_id, genome_pos);
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

    cur_genome_threads_[genome_id]->push_back(cur_edge);
    thread.push_back(make_pair(graph_pos, genome_pos));
    graph_pos += g_->length(cur_edge);
    genome_pos = edge_ranges_[cur_edge].GetForwardPos(cur_pos);
    cur_pos.second = genome_pos;

    const VertexId v = g_->EdgeEnd(cur_edge);

    DEBUG("current edge " << g_->str(cur_edge) << ", outgoing count " << g_->OutgoingEdgeCount(v));
    cur_edge = EdgeId(0);
    for (const auto &out_edge : g_->OutgoingEdges(v)) {
      DEBUG("considering edge " << g_->str(out_edge) << " at position (seq) " << genome_pos);

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
CoordinatesHandler<Graph>::FindGenomeFirstEdge(const uchar genome_id) const {
  std::pair<uchar, size_t> pos_in_question(genome_id, 0);

  for (auto it = g_->SmartEdgeBegin(); !it.IsEnd(); ++it) {
    auto range_it = edge_ranges_.find(*it);
    if (range_it == edge_ranges_.end())
      continue;
    if (range_it->second.HasForwardLink(pos_in_question)) {
      return *it;
    }
  }

  // remember first edge and update it
  throw std::runtime_error("Could not find start of the sequence in graph");
  return EdgeId(0);
}

template <class Graph>
size_t CoordinatesHandler<Graph>::GetOriginalPos(
    const uchar genome_id, const size_t new_pos) const {
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

    // Kmers can have different lengths so going from larger kmers to smaller
    // implies shorting of thread length what may lead to "range-overflow"
    if (cur_pos > thread_it->back().first)
      cur_pos = thread_it->back().first;

    auto found_it = std::lower_bound(thread_it->begin(), thread_it->end(),
                                     make_pair(cur_pos, size_t(0)));

    DEBUG("Searching for pos " << cur_pos << "in thread of " << thread_it->front() << " - " << thread_it->back());
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

    DEBUG("from ranges " << graph_range << " and " << genome_range << " in search of " << cur_pos);
    cur_pos = CalculatePos(graph_range, genome_range, cur_pos);
  }

  return cur_pos;
}

template<class Graph>
size_t CoordinatesHandler<Graph>::GetNewestPos(
    const uchar genome_id, const size_t old_pos) const {
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

  DEBUG("from ranges " << graph_range << " and " << genome_range << " in search of " << search_pos);
  return CalculatePos(graph_range, genome_range, search_pos);
}

/*
 *  Some handling methods (description in omni_utils.hpp)
 */
template <class Graph>
void CoordinatesHandler<Graph>::HandleDelete(EdgeId e) {
  /*
  if (HasEdgeData(e)) {
    INFO("edge " << g_->str(e) << " " << edge_ranges_[e].DebugOutput());
  }
  */
  VERIFY(!HasEdgeData(e));
  CleanEdgeData(e);
}

template <class Graph>
void CoordinatesHandler<Graph>::HandleMerge(const vector<EdgeId> &old_edges, EdgeId new_edge) {
  //TRACE("HandleMerge : " << Debug(old_edges) << " -> " << g_->str(new_edge));
  for (const auto &edge : old_edges) {
    if (HasEdgeData(edge)) {
      edge_ranges_[new_edge] += edge_ranges_[edge];
      CleanEdgeData(edge);
    }
  }
}

template <class Graph>
void CoordinatesHandler<Graph>::HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
  //TRACE("HandleGlue : " << g_->str(new_edge) << " <- " << g_->str(edge1) << " + " << g_->str(edge2));
  if (HasEdgeData(edge1)) {
    edge_ranges_[new_edge] += edge_ranges_[edge1];
    CleanEdgeData(edge1);
  }
  if (HasEdgeData(edge2)) {
    edge_ranges_[new_edge] += edge_ranges_[edge2];
    CleanEdgeData(edge2);
  }
}

template <class Graph>
void CoordinatesHandler<Graph>::HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
                         EdgeId new_edge_2) {
  //TRACE("HandleSplit " << old_edge << " -> " << new_edge_1 << " + " << new_edge_2);
  VERIFY(!HasEdgeData(old_edge));
  /*
  const std::vector<std::pair<uchar, Range> > old_ranges =
      edge_ranges[old_edge].GetRanges();
  for (const auto &range : old_ranges) {

  }
  */
}

}
