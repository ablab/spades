//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "xmath.h"
#include "omni/omni_utils.hpp"
#include "sequence/sequence.hpp"

#include <boost/iterator/iterator_facade.hpp>

#include <btree/btree_set.h>

#include <cmath>
#include <map>
#include <limits>

//#define MERGE_DATA_RELATIVE_DIFFERENCE 0.3

namespace omnigraph {

namespace de {

/**
 * PairInfo class represents basic data unit for paired information: edges first and second appear
 * in genome at distance d_ and this information has weight weight_.
 */
struct Point {
  public:
    typedef double value_type;

    value_type d;
    value_type weight;
    value_type var;

    Point()
            : d(0.0), weight(0.0), var(0.0) {}

    explicit Point(value_type distance, value_type weight, value_type variance)
            : d(distance), weight(weight), var(variance) {}

    Point(const Point& rhs)
            : d(rhs.d), weight(rhs.weight), var(rhs.var) {}

    std::string str() const {
        stringstream ss;
        ss << "Point: " << " distance = " << this->d
           << ", weight = " << this->weight
           << ", variance = " << this->var;
        return ss.str();
    }

    Point& operator=(const Point& rhs) {
        using namespace math;
        update_value_if_needed<value_type>(d, rhs.d);
        update_value_if_needed<value_type>(weight, rhs.weight);
        update_value_if_needed<value_type>(var, rhs.var);
        return *this;
    }

    bool operator<(const Point& rhs) const {
        return math::ls(this->d, rhs.d);
    }

    bool operator==(const Point& rhs) const {
        return math::eq(this->d, rhs.d);
    }

    bool operator!=(const Point& rhs) const {
        return !(operator==(rhs));
    }

    Point operator-() const {
        return Point(-d, weight, var);
    }

    Point operator+(const Point &rhs) const {
        value_type weight_rhs = rhs.weight;
        // counting new bounds in the case, when we are merging pair infos with var != 0
        value_type left_bound = std::min(d - var, rhs.d - rhs.var);
        value_type right_bound = std::max(d + var, rhs.d + rhs.var);
        value_type new_dist = (left_bound + right_bound) * 0.5;
        value_type new_weight = weight + weight_rhs;
        value_type new_variance = (right_bound - left_bound) * 0.5;

        return Point(new_dist, new_weight, new_variance);
    }
};

inline int rounded_d(Point p) {
    return math::round_to_zero(p.d);
}

inline std::ostream& operator<<(std::ostream& os, const Point &point) {
    return os << point.str();
}

//typedef std::set<Point> Histogram;
typedef btree::btree_set<Point> Histogram;

inline bool ClustersIntersect(Point p1, Point p2) {
  return math::le(p1.d, p2.d + p1.var + p2.var) &&
         math::le(p2.d, p1.d + p1.var + p2.var);
}

inline Point ConjugatePoint(size_t l1, size_t l2, const Point& point) {
  return Point(point.d + (double) l2 - (double) l1, point.weight, point.var);
}

// tuple of a pair of edges @first, @second, and a @point
template<typename EdgeId>
struct PairInfo {
    EdgeId first;
    EdgeId second;
    Point point;

    PairInfo()
            : first(), second(), point() {}


    PairInfo(const PairInfo& pair_info)
            : first(pair_info.first), second(pair_info.second), point(pair_info.point)
    {}

    PairInfo(EdgeId first, EdgeId second, double d, double weight, double var)
            : first(first), second(second), point(d, weight, var)
    {}

    PairInfo(EdgeId first, EdgeId second, Point point)
            : first(first), second(second), point(point) {}

    // Two paired infos are considered equal
    // if they coinside in all parameters except for weight and variance.
    bool operator==(const PairInfo& rhs) const {
        const PairInfo &lhs = *this;
        return lhs.first == rhs.first && lhs.second == rhs.second && lhs.point == rhs.point;
    }

    bool operator!=(const PairInfo& rhs) const {
        return !(*this == rhs);
    }

    bool operator<(const PairInfo<EdgeId>& rhs) const {
        const PairInfo<EdgeId>& lhs = *this;
        return (lhs.first == rhs.first ?
                (lhs.second == rhs.second ? lhs.point < rhs.point : lhs.second < rhs.second)
                : lhs.first  < rhs.first);
    }

    double d() const      { return point.d;      }
    double weight() const { return point.weight; }
    double var() const    { return point.var;    }

    double& d()           { return point.d;      }
    double& weight()      { return point.weight; }
    double& var()         { return point.var;    }
};

template<typename EdgeId>
ostream& operator<<(ostream& os, const PairInfo<EdgeId>& info) {
  return os << "PairInfo: first = " << info.first << ", second = " << info.second
            << "Point : " << info.point;
}

template<typename EdgeId>
const PairInfo<EdgeId> MinPairInfo(EdgeId id) {
  return PairInfo<EdgeId>(id, EdgeId(typename EdgeId::pointer_type(1)),
      -10000000000, 0., 0.);
}

template<typename EdgeId>
const PairInfo<EdgeId> MaxPairInfo(EdgeId id) {
  return PairInfo<EdgeId>(id, EdgeId(typename EdgeId::pointer_type(-1)),
       10000000000, 0., 0.);
}

template<typename EdgeId>
const PairInfo<EdgeId> MinPairInfo(EdgeId e1, EdgeId e2) {
  PairInfo<EdgeId> info = MinPairInfo(e1);
  info.second = e2;
  return info;
}

template<typename EdgeId>
const PairInfo<EdgeId> MaxPairInfo(EdgeId e1, EdgeId e2) {
  PairInfo<EdgeId> info = MaxPairInfo(e1);
  info.second = e2;
  return info;
}

/**
 * Method returns approximate distance between occurrences of edges in genome rounded to the nearest
 * integer. In case of a tie closest to 0 value is chosen thus one can assume that distance
 * is rounded the same way as opposite one.
 * todo check that written here is true
 */
template<typename EdgeId>
inline int rounded_d(PairInfo<EdgeId> const& pi) {
  return math::round_to_zero(pi.d());
}

template<typename EdgeId>
PairInfo<EdgeId> BackwardInfo(const PairInfo<EdgeId>& pi) {
  return PairInfo<EdgeId>(pi.second, pi.first, -pi.point);
}

template<typename EdgeId>
inline bool IsSymmetric(PairInfo<EdgeId> const& pi) {
  return pi.first == pi.second && math::eq(pi.d(), 0.);
}

// new map { EdgeId -> (EdgeId -> (d, weight, var)) }
template<class Graph,
         class InnerMapType = std::map<typename Graph::EdgeId, Histogram>,
         class IndexDataType = std::map<typename Graph::EdgeId, InnerMapType> >
class PairedInfoStorage {
 public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Histogram::const_iterator HistIterator;
    typedef typename Point::value_type PointValueType;
    typedef InnerMapType InnerMap;
    typedef typename IndexDataType::const_iterator DataIterator;

    PairedInfoStorage()
            : size_(0) {}

    class EdgePairIterator :
            public boost::iterator_facade<EdgePairIterator,
                                          const Histogram,
                                          boost::forward_traversal_tag,
                                          const Histogram& > {

      public:
        EdgePairIterator(DataIterator cedge, DataIterator eedge)
                : cedge_(cedge), eedge_(eedge), sedge_() {
            if (cedge_ == eedge_)
                return;

            sedge_ = cedge_->second.begin();
            skip_empty();
        }

        EdgeId first() const { return cedge_->first; }
        EdgeId second() const { return sedge_->first; }

        friend ostream& operator<<(ostream& os, const EdgePairIterator& iter) {
            return os << iter.first() << " " << iter.second();
        }

      private:
        typedef typename InnerMap::const_iterator InnerIterator;

        friend class boost::iterator_core_access;

        void skip_empty() {
            while (sedge_ == cedge_->second.end()) {
                ++cedge_;
                if (cedge_ == eedge_)
                    break;
                sedge_ = cedge_->second.begin();
            }
        }

        void increment() {
            ++sedge_;
            skip_empty();
        }

        bool equal(const EdgePairIterator &other) const {
            return other.cedge_ == cedge_ && (cedge_ == eedge_ || other.sedge_ == sedge_);
        }

        const Histogram& dereference() const {
            return sedge_->second;
        }

        DataIterator cedge_, eedge_;
        InnerIterator sedge_;
    };

    class EdgeIterator :
            public boost::iterator_facade<EdgeIterator,
                                          const std::pair<EdgeId, Point>,
                                          boost::forward_traversal_tag,
                                          const std::pair<EdgeId, Point> > {
        typedef typename Histogram::const_iterator histogram_iterator;
        typedef typename InnerMap::const_iterator InnerIterator;

      public:
        EdgeIterator(InnerIterator cedge, InnerIterator eedge)
                : cedge_(cedge), eedge_(eedge), point_() {
            if (cedge_ == eedge_)
                return;

            point_ = cedge_->second.begin();
            skip_empty();
        }

      private:
        friend class boost::iterator_core_access;

        void skip_empty() {
            while (point_ == cedge_->second.end()) {
                ++cedge_;
                if (cedge_ == eedge_)
                    break;
                point_ = cedge_->second.begin();
            }
        }

        void increment() {
            ++point_;
            skip_empty();
        }

        bool equal(const EdgeIterator &other) const {
            return other.cedge_ == cedge_ && (cedge_ == eedge_ || other.point_ == point_);
        }

        const std::pair<EdgeId, Point> dereference() const {
            return std::make_pair(cedge_->first, *point_);
        }

        InnerIterator cedge_, eedge_;
        histogram_iterator point_;
    };

    EdgePairIterator begin() const {
        return EdgePairIterator(index_.begin(), index_.end());
    }

    EdgePairIterator end() const {
        return EdgePairIterator(index_.end(), index_.end());
    }

    DataIterator Begin() const {
        return index_.begin();
    }

    DataIterator End() const {
        return index_.end();
    }

    EdgeIterator edge_begin(EdgeId edge) const {
        VERIFY(contains(edge));
        return edge_begin(index_.find(edge));
    }

    EdgeIterator edge_end(EdgeId edge) const {
        VERIFY(contains(edge));
        return edge_end(index_.find(edge));
    }

    bool contains(EdgeId edge) const {
        return index_.count(edge);
    }

    // FIXME: Make these private
    EdgeIterator edge_begin(typename IndexDataType::const_iterator entry) const {
        return EdgeIterator(entry->second.begin(), entry->second.end());
    }

    EdgeIterator edge_end(typename IndexDataType::const_iterator entry) const {
        return EdgeIterator(entry->second.end(), entry->second.end());
    }

    // adding pair infos
    void AddPairInfo(const std::pair<EdgeId, EdgeId>& edge_pair,
                     Point point_to_add,
                     bool add_reversed = true) {
        AddPairInfo(edge_pair.first, edge_pair.second, point_to_add, add_reversed);
    }

    void AddPairInfo(const std::pair<EdgeId, EdgeId>& edge_pair,
                     PointValueType d, PointValueType weight, PointValueType var,
                     bool add_reversed = true) {
        AddPairInfo(edge_pair.first, edge_pair.second, Point(d, weight, var), add_reversed);
    }

    void AddPairInfo(EdgeId e1, EdgeId e2,
                     PointValueType d, PointValueType weight, PointValueType var,
                     bool add_reversed = true) {
        AddPairInfo(e1, e2, Point(d, weight, var), add_reversed);
    }

    void AddPairInfo(EdgeId e1, EdgeId e2,
                     Point point_to_add,
                     bool add_reversed = true) {
        Histogram& histogram = index_[e1][e2];
        HistIterator iterator_to_point = histogram.find(point_to_add);

        if (iterator_to_point != histogram.end())
            MergeData(e1, e2, *iterator_to_point, point_to_add, add_reversed);
        else
            InsertPoint(e1, e2, histogram, point_to_add, add_reversed);
    }

    void DeletePairInfo(EdgeId e1, EdgeId e2,
                        Point point_to_remove) {
        Histogram& histogram = index_[e1][e2];
        histogram.erase(point_to_remove);
    }

    // erasing specific entry from the index
    size_t RemovePairInfo(EdgeId e1, EdgeId e2, const Point& point_to_remove) {
        auto iter = index_.find(e1);
        if (iter != index_.end()) {
            InnerMap& map = iter->second;
            auto iter2 = map.find(e2);
            if (iter2 != map.end()) {
                Histogram& hist = iter2->second;
                size_t success = hist.erase(point_to_remove);
                if (success == 1)
                    --size_;
                if (hist.empty())
                    map.erase(iter2);
                if (map.empty())
                    index_.erase(iter);

                return success;
            }
        }
        return 0;
    }

    void RemovePairInfo(const PairInfo<EdgeId>& info) {
        this->RemovePairInfo(info.first, info.second, info.point);
    }

    // removing all points from @e1, @e2 histogram,
    // returns 1 if operation was successful, 0 if not
    size_t RemoveEdgePairInfo(EdgeId e1, EdgeId e2) {
        auto iter = index_.find(e1);
        if (iter != index_.end()) {
            InnerMap& map = iter->second;
            auto iter2 = map.find(e2);
            if (iter2 != map.end()) {
                Histogram& hist = iter2->second;
                size_t size_decrease = hist.size();
                map.erase(iter2);
                size_ -= size_decrease;
                if (map.empty())
                    index_.erase(iter);

                return 1;
            }
        }
        return 0;
    }

    // removes all points, which refer to this edge
    // also removes all backward information
    void RemoveEdgeInfo(EdgeId edge) {
        InnerMap& inner_map = index_[edge];
        for (auto iter = inner_map.begin(); iter != inner_map.end(); ++iter) {
            EdgeId e2 = iter->first;
            if (edge != e2)
                this->RemoveEdgePairInfo(e2, edge);
        }
        size_t size_of_removed = inner_map.size();
        index_.erase(edge);
        size_ -= size_of_removed;
    }

    void Clear() {
        index_.clear();
        size_ = 0;
    }

    // Usual implementation, the same as in the old paired index
    std::vector<PairInfo<EdgeId> > GetEdgeInfo(EdgeId edge) const {
        typename IndexDataType::const_iterator iter = index_.find(edge);
        TRACE("Getting edge info");
        if (iter == index_.end())
            return std::vector<PairInfo<EdgeId> >();

        std::vector<PairInfo<EdgeId> > result;
        result.reserve(iter->second.size());
        for (auto I = edge_begin(iter), E = edge_end(iter); I != E; ++I) {
            std::pair<EdgeId, Point> entry = *I;
            result.push_back(PairInfo<EdgeId>(edge, entry.first, entry.second));
        }
        return result;
    }

    // faster implementation, but less resolver-friendly
    // returns InnerMap instead of vector<>,
    // one can iterate it using FastIterator class
    const InnerMap GetEdgeInfo(EdgeId edge, int) const {
        typename IndexDataType::const_iterator iter = index_.find(edge);
        if (iter == index_.end())
            return InnerMap();
        else
            return iter->second;
    }

    const Histogram GetEdgePairInfo(EdgeId e1, EdgeId e2) const {
        typename IndexDataType::const_iterator iter = index_.find(e1);
        if (iter == index_.end())
            return Histogram();
        else {
            const InnerMap& inner_map = iter->second;
            typename InnerMap::const_iterator iter2 = inner_map.find(e2);
            if (iter2 == inner_map.end())
                return Histogram();
            else
                return iter2->second;
        }
    }

    void Prune() {
        for (auto iter = index_.begin(); iter != index_.end(); ) {
            // First, remove all the empty Histograms
            InnerMap& inner_map = iter->second;
            for (auto it = inner_map.begin(); it != inner_map.end(); ) {
                if (it->second.empty())
                    inner_map.erase(it++);
                else
                    ++it;
            }

            // Now, pretty much the same, but the outer stuff
            if (inner_map.empty())
                index_.erase(iter++);
            else
                ++iter;
        }
    }

    // here we trying to insert PairInfo,
    // if there is no existing PairInfo with the same key
    // very complicated implementation, but it seems to be faster.
    template<class Storage>
    void AddAll(const Storage& index_to_add) {
        IndexDataType& base_index = this->index_;
        for (auto AddI = index_to_add.Begin(), E = index_to_add.End(); AddI != E; ++AddI) {
            EdgeId e1_to_add = AddI->first;
            const auto& map_to_add = AddI->second;
            InnerMap& map_already_exists = base_index[e1_to_add];
            MergeInnerMaps(e1_to_add, map_to_add, map_already_exists);
        }
    }

    size_t size() const {
        return size_;
    }

  private:
    bool IsSymmetric(EdgeId e1, EdgeId e2,
                     Point point) const {
        return (e1 == e2) && math::eq(point.d, 0.);
    }

    // modifying the histogram
    void InsertPoint(EdgeId e1, EdgeId e2,
                     Histogram& histogram,
                     Point new_point,
                     bool add_reversed) {
        // first backwards
        if (add_reversed && !IsSymmetric(e1, e2, new_point)) {
            index_[e2][e1].insert(-new_point);
            ++size_;
        }

        histogram.insert(new_point);
        ++size_;
    }

    void UpdateSinglePoint(Histogram &hist, Histogram::iterator point_to_update, Point new_point) {
        Histogram::iterator after_removed = hist.erase(point_to_update);
        hist.insert(after_removed, new_point);
    }

    void MergeData(EdgeId e1, EdgeId e2,
                   Point point_to_update,
                   Point point_to_add,
                   bool add_reversed) {
        if (add_reversed) {
            Histogram& histogram = index_[e2][e1];
            UpdateSinglePoint(histogram, histogram.find(-point_to_update), -(point_to_update + point_to_add));
        }

        Histogram& histogram = index_[e1][e2];
        UpdateSinglePoint(histogram, histogram.find(point_to_update), point_to_update + point_to_add);
    }

    void MergeData(Histogram& hist,
                   Histogram::iterator to_update,
                   Point point_to_add) {
        UpdateSinglePoint(hist, to_update, *to_update + point_to_add);
    }

    template<class OtherMap>
    void MergeInnerMaps(EdgeId /*e1_to_add*/,
                        const OtherMap& map_to_add,
                        InnerMap& map) {
        typedef typename Histogram::iterator hist_iterator;
        typedef typename InnerMap::iterator map_iterator;
        for (auto I = map_to_add.begin(), E = map_to_add.end(); I != E; ++I) {
            const Histogram& hist_to_add = I->second;
            const pair<map_iterator, bool>& result = map.insert(*I);
            if (!result.second) { // in this case we need to merge two hists
                Histogram& hist_exists = (result.first)->second;
                // pretty much the same
                for (auto p_it = hist_to_add.begin(), E = hist_to_add.end(); p_it != E; ++p_it) {
                    Point new_point = *p_it;
                    const pair<hist_iterator, bool>& result = hist_exists.insert(new_point);
                    if (!result.second) { // in this case we need to merge two points
                        MergeData(hist_exists, result.first, new_point);
                    } else
                        ++size_;
                }
            } else
                size_ += hist_to_add.size();
        }
    }

  protected:

    void TransferInfo(EdgeId old_edge, EdgeId new_edge,
                      int shift = 0,
                      double weight_scale = 1.) {
        const InnerMap& inner_map = this->GetEdgeInfo(old_edge, 0);
        for (auto iter = inner_map.begin(); iter != inner_map.end(); ++iter) {
            EdgeId e2 = iter->first;
            const Histogram& histogram = iter->second;
            for (auto point_iter = histogram.begin(); point_iter != histogram.end(); ++point_iter) {
                Point cur_point = *point_iter;
                if (old_edge != e2) {
                    AddPairInfo(new_edge, e2,
                                cur_point.d - shift,
                                weight_scale * cur_point.weight,
                                cur_point.var);
                } else if (!math::eq(cur_point.d, 0.)) {
                    AddPairInfo(new_edge, new_edge,
                                cur_point.d,
                                weight_scale * 0.5 * cur_point.weight,
                                cur_point.var);
                } else {
                    AddPairInfo(new_edge, new_edge,
                                cur_point.d,
                                weight_scale * cur_point.weight,
                                cur_point.var);
                }
            }
        }
    }

  private:
    IndexDataType index_;
    size_t size_;
};

template<class Graph>
using PairedInfoBuffer = PairedInfoStorage<Graph,
                                           std::unordered_map<typename Graph::EdgeId, Histogram>,
                                           std::unordered_map<typename Graph::EdgeId,
                                                              std::unordered_map<typename Graph::EdgeId, Histogram> > >;

template<class Graph>
class PairedInfoIndexT: public GraphActionHandler<Graph>, public PairedInfoStorage<Graph> {
  public:
    typedef typename Graph::EdgeId EdgeId;

    PairedInfoIndexT(const Graph& graph) :
            GraphActionHandler<Graph>(graph, "PairedInfoIndexT") {}

    ~PairedInfoIndexT() {
        TRACE("~PairedInfoIndexT ok");
    }

    void Init() {
        for (auto it = this->g().ConstEdgeBegin(); !it.IsEnd(); ++it) {
            this->HandleAdd(*it);
        }
    }

    // method adds paired info to the conjugate edges
    void AddConjPairInfo(EdgeId e1, EdgeId e2,
                         Point point_to_add,
                         bool add_reversed = 1) {
        const Graph& g = this->g();
        this->AddPairInfo(g.conjugate(e2),
                          g.conjugate(e1),
                          ConjugatePoint(g.length(e1), g.length(e2), point_to_add),
                          add_reversed);
    }

    // prints the contents of index
    void PrintAll() const {
        size_t size = 0;
        for (auto I = this->begin(), E = this->end(); I != E; ++I) {
            EdgeId e1 = I.first();
            EdgeId e2 = I.second();
            const Histogram& histogram = *I;
            size += histogram.size();
            INFO("Histogram for edges "
                 << this->g().int_id(e1) << " "
                 << this->g().int_id(e2));
            for (auto it = histogram.begin(); it != histogram.end(); ++it) {
                INFO("    Entry " << it->str());
            }
        }
        VERIFY_MSG(this->size() == size, "Size " << size << " must have been equal to " << this->size());
    }

    // Handlers
    virtual void HandleAdd(EdgeId edge) {
        const Graph& graph = this->g();
#pragma omp critical
        {
            TRACE("Handling Addition " << graph.int_id(edge));
            this->AddPairInfo(edge, edge, 0., 0., 0.);
        }
    }

    virtual void HandleDelete(EdgeId edge) {
        const Graph& graph = this->g();

        TRACE("Handling Deleting " << graph.int_id(edge));
        this->RemoveEdgeInfo(edge);
    }

    virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
        TRACE("Handling Merging");
        this->AddPairInfo(new_edge, new_edge, 0., 0., 0.);
        int shift = 0;
        const Graph& graph = this->g();
        for (size_t i = 0; i < old_edges.size(); ++i) {
            EdgeId old_edge = old_edges[i];
            this->TransferInfo(old_edge, new_edge, shift);
            shift -= (int) graph.length(old_edge);
        }
    }

    virtual void HandleGlue(EdgeId new_edge, EdgeId e1, EdgeId e2) {
        const Graph& graph = this->g();
        TRACE("Handling Glueing " << graph.int_id(new_edge) << " " << graph.int_id(e1) << " "
              << graph.int_id(e2));
        this->TransferInfo(e2, new_edge);
        this->TransferInfo(e1, new_edge);
    }

    virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge_2) {
        const Graph& graph = this->g();
        TRACE("Handling Splitting " << graph.int_id(old_edge) << " " << graph.int_id(new_edge_1)
              << " " << graph.int_id(new_edge_2));
        double ratio = (double) graph.length(new_edge_1) * 1. / (double) graph.length(old_edge);
        this->TransferInfo(old_edge, new_edge_1, 0, ratio);
        this->TransferInfo(old_edge, new_edge_2, (int) graph.length(new_edge_1), 1. - ratio);
    }

  DECL_LOGGER("PairedInfoIndexT");
};

template <class Graph>
struct PairedInfoIndicesT {
    typedef PairedInfoIndexT<Graph> IndexT;
    std::vector<IndexT> data_;

    PairedInfoIndicesT(const Graph& graph, size_t lib_num) {
        for (size_t i = 0; i < lib_num; ++i)
            data_.emplace_back(graph);
    }

    void Init() {
        for (auto& it : data_)
            it.Init();
    }

    void Attach() {
        for (auto& it : data_)
            it.Attach();
    }

    void Detach() {
        for (auto& it : data_)
            it.Detach();
    }

    IndexT& operator[](size_t i) {
        return data_[i];
    }

    const IndexT& operator[](size_t i) const {
        return data_[i];
    }

    size_t size() const {
        return data_.size();
    }
};

//New metric weight normalizer
template<class Graph>
class PairedInfoWeightNormalizer {
  typedef typename Graph::EdgeId EdgeId;
  const Graph& g_;
  const size_t insert_size_;
  //todo use this param!
  const double is_var_;
  const size_t read_length_;
  const size_t k_;
  const double avg_coverage_;
public:

  //Delta better to be around 5-10% of insert size
  PairedInfoWeightNormalizer(const Graph& g, size_t insert_size,
      double is_var, size_t read_length, size_t k, double avg_coverage) :
      g_(g), insert_size_(insert_size), is_var_(is_var), read_length_(
          read_length), k_(k), avg_coverage_(avg_coverage) {
  }

  const PairInfo<EdgeId> NormalizeWeightWithCoverage(const PairInfo<EdgeId>& pair_info) {
      PairInfo<EdgeId> new_info = pair_info;
      new_info.weight() *= g_.length(pair_info.first) * g_.length(pair_info.second) * 1.
                        / (g_.coverage(pair_info.first) * g_.coverage(pair_info.second));
      return new_info;
  }

  const Point NormalizeWeight(EdgeId e1, EdgeId e2, Point point) const {
    double w = 0.;
    if (math::eq(point.d, 0.) && e1 == e2) {
      w = 0. + (double) g_.length(e1) - (double) insert_size_ + 2. * (double) read_length_ + 1. - (double) k_;
    }
    else {
      if (math::ls(point.d, 0.)) {
        using std::swap;
        swap(e1, e2);
      }
      int gap_len = abs(rounded_d(point)) - (int) g_.length(e1);
      int right = std::min((int) insert_size_, gap_len + (int) g_.length(e2) + (int) read_length_);
      int left = std::max(gap_len, (int) insert_size_ - (int) read_length_ - (int) g_.length(e1));
      w = 0. + (double) (right - left + 1 - (int) k_);
    }

    double result_weight = point.weight;
    if (math::gr(w, /*-10.*/0.)) {
      result_weight /= w; //(w + 10);
    } else
      result_weight = 0.;

    double cov_norm_coeff = avg_coverage_ / (2. * (double) (read_length_ - k_));
    result_weight /= cov_norm_coeff;

    Point result(point);
    result.weight = result_weight;
    return result;
  }
};

template<class Graph>
class JumpingNormalizerFunction {
private:
  typedef typename Graph::EdgeId EdgeId;
  const Graph& graph_;
  size_t read_length_;
  size_t max_norm_;

public:
  JumpingNormalizerFunction(const Graph& graph, size_t read_length,
      size_t max_norm) :
      graph_(graph), read_length_(read_length), max_norm_(max_norm) {
  }

  size_t norm(EdgeId e1, EdgeId e2) const {
    return std::min(std::min(graph_.length(e1), graph_.length(e2)),
        max_norm_) + read_length_ - graph_.k();
  }

  const Point operator()(EdgeId e1, EdgeId e2, Point point) const {
    return Point(point.d,
                 point.weight * 1. / (double) norm(e1, e2),
                 point.var);
  }
};

template<class Graph>
const Point TrivialWeightNormalization(typename Graph::EdgeId,
                                       typename Graph::EdgeId,
                                       Point point)
{
  return point;
}

template<class Graph>
class PairedInfoNormalizer {
  typedef typename Graph::EdgeId EdgeId;

 public:
  typedef boost::function<const Point(EdgeId, EdgeId, Point)> WeightNormalizer;

  PairedInfoNormalizer(WeightNormalizer normalizing_function) :
      normalizing_function_(normalizing_function)
  {
  }

// temporary due to path_extend absolute thresholds
  void FillNormalizedIndex(const PairedInfoIndexT<Graph>& paired_index,
                                 PairedInfoIndexT<Graph>& normalized_index,
                                 double coeff = 1.) const
  {
    for (auto I = paired_index.begin(), E = paired_index.end(); I != E; ++I) {
      const Histogram& hist = *I;
      EdgeId e1 = I.first();
      EdgeId e2 = I.second();
      TRACE("first second " << e1 << " " << e2);
      for (auto it2 = hist.begin(); it2 != hist.end(); ++it2) {
        Point tmp(*it2);
        TRACE("TEMP point " << tmp);
        tmp.weight *= coeff;
        TRACE("Normalized pair info " << tmp << " " << normalizing_function_(e1, e2, tmp));
        normalized_index.AddPairInfo(e1, e2, normalizing_function_(e1, e2, tmp), false);
      }
    }
  }

 private:
  WeightNormalizer normalizing_function_;
  DECL_LOGGER("PairedInfoNormalizer");
};

};

}
