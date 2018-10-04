//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * edges_position_handler.hpp
 *
 *  Created on: 22.07.2011
 *
 */

#ifndef EDGES_POSITION_HANDLER_HPP_
#define EDGES_POSITION_HANDLER_HPP_

#include "utils/stl_utils.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/core/action_handlers.hpp"
#include <map>

namespace omnigraph {

struct EdgePosition {
    std::string contigId;
    MappingRange mr;
    EdgePosition(const std::string &_contigId, MappingRange _mr) : contigId(_contigId), mr(_mr) {
    }

    EdgePosition() {
    }
};

inline std::ostream &operator<<(std::ostream &os, const EdgePosition &ep) {
    return os << ep.contigId << " " << ep.mr;
}

template<class Graph>
class EdgesPositionHandler: public GraphActionHandler<Graph> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef std::set<MappingRange> RangeSet;

    size_t max_mapping_gap_;
    size_t max_gap_diff_;
    std::map<EdgeId, std::map<std::string, RangeSet>> edges_positions_;
    //TODO extract RangeSet aka set<MappingRange> as a storage class

    MappingRange EraseAndExtract(RangeSet &ranges, RangeSet::iterator position, const MappingRange &new_pos) const {
        auto old_pos = *position;
        if (old_pos.IntersectLeftOf(new_pos) || old_pos.StrictlyContinuesWith(new_pos, max_mapping_gap_, max_gap_diff_)) {
            ranges.erase(position);
            return old_pos.Merge(new_pos);
        } else if (new_pos.IntersectLeftOf(old_pos) || new_pos.StrictlyContinuesWith(old_pos, max_mapping_gap_, max_gap_diff_)) {
            ranges.erase(position);
            return new_pos.Merge(old_pos);
        } else {
            return new_pos;
        }
    }

    std::string RangeStr(const Range &range) const {
        std::stringstream ss;
        ss << "[" << (range.start_pos + 1) << " - " << range.end_pos << "]";
        return ss.str();
    }

public:
    MappingRange EraseAndExtract(RangeSet &ranges, MappingRange new_pos) const {
        auto it = ranges.lower_bound(new_pos);
        if (it != ranges.end()) {
            new_pos = EraseAndExtract(ranges, it, new_pos);
            it = ranges.lower_bound(new_pos);
        }
        if (it != ranges.begin()) {
            new_pos = EraseAndExtract(ranges, std::prev(it), new_pos);
        }
        return new_pos;
    }

    RangeSet GetEdgePositions(EdgeId edge, const std::string &contig_id) const {
        VERIFY(this->IsAttached());
        auto edge_it = edges_positions_.find(edge);
        if (edge_it == edges_positions_.end())
            return {};
        const auto &positions = edge_it->second;
        auto it = positions.find(contig_id);
        if (it == positions.end())
            return {};
        else
            return it->second;
    }

    MappingRange GetUniqueEdgePosition(EdgeId edge, const std::string &contig_id) const {
        auto poss = GetEdgePositions(edge, contig_id);
        VERIFY(poss.size() == 1);
        return *poss.begin();
    }

    std::vector<EdgePosition> GetEdgePositions(EdgeId edge) const {
        VERIFY(this->IsAttached());
        auto edge_it = edges_positions_.find(edge);
        if (edge_it == edges_positions_.end())
            return {};
        std::vector<EdgePosition> result;
        for (const auto &i : edge_it->second) {
            for (const auto &pos : i.second) {
                result.push_back(EdgePosition(i.first, pos));
            }
        }
        return result;
    }

    void AddEdgePosition(EdgeId edge, const std::string &contig_id, size_t start, size_t end, size_t m_start, size_t m_end) {
        VERIFY(this->IsAttached());
        AddEdgePosition(edge, contig_id, MappingRange(start, end, m_start, m_end));
    }

    void AddEdgePosition(EdgeId edge, const std::string &contig_id, MappingRange new_pos) {
        VERIFY(this->IsAttached());
        if (new_pos.empty())
            return;
        auto &new_set = edges_positions_[edge][contig_id];
        new_pos = EraseAndExtract(new_set, new_pos);
        new_set.insert(new_pos);
    }

    void AddAndShiftEdgePositions(EdgeId edge, const std::map<std::string, RangeSet> &contig_map, int shift = 0) {
        VERIFY(this->IsAttached());
        for (const auto &contig : contig_map) {
            for (const auto &e : contig.second) {
                AddEdgePosition(edge, contig.first, e.Shift(shift).Fit(this->g().length(edge)));
            }
        }
    }

    template<typename Iter>
    void AddEdgePositions(EdgeId edge, Iter begin, Iter end) {
        VERIFY(this->IsAttached());
        for (auto it = begin; it != end; ++it) {
            AddEdgePosition(edge, it->contigId, it->mr);
        }
    }

    std::string str(EdgeId edge) const {
        VERIFY(this->IsAttached());
        std::vector<EdgePosition> positions = GetEdgePositions(edge);
        size_t counter = 0;
        std::stringstream ss;
        for (const auto &pos : positions) {
            ss << "(" << pos.contigId << ": "
               << RangeStr(pos.mr.initial_range) << " --> "
               << RangeStr(pos.mr.mapped_range) << ")\\n";
            counter++;
            if (counter > 30) {
                ss << "and many more. Totally " << positions.size() << " positions.";
                break;
            }
        }
        return ss.str();
    }
    
    /**
    * @param max_mapping_gap - maximal difference in positions of 
    * original sequence for two mapping ranges to be merged.
    * @param max_gap_diff - maximal difference between gaps in initial and mapped ranges for
    * mapping ranges to be merged
    */
    EdgesPositionHandler(const Graph &g, size_t max_mapping_gap, size_t max_gap_diff = 0) :
            GraphActionHandler<Graph>(g, "EdgePositionHandler"), 
            max_mapping_gap_(max_mapping_gap),
            max_gap_diff_(max_gap_diff) {
    }

    virtual ~EdgesPositionHandler() {
        TRACE("~EdgePositionHandler ok");
    }

    virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
//        TRACE("Handle glue ");
        auto positions1 = GetEdgePositions(edge1);
        auto positions2 = GetEdgePositions(edge2);
        AddEdgePositions(new_edge, positions1.begin(), positions1.end());
        AddEdgePositions(new_edge, positions2.begin(), positions2.end());
    }

    virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
        if (oldEdge == this->g().conjugate(oldEdge)) {
            WARN("EdgesPositionHandler does not support self-conjugate splits");
            return;
        }
        if (edges_positions_.count(oldEdge) != 0) {
            auto contig_map = edges_positions_[oldEdge];
            AddAndShiftEdgePositions(newEdge1, contig_map, 0);
            AddAndShiftEdgePositions(newEdge2, contig_map, -int(this->g().length(newEdge1)));
        }
    }

    virtual void HandleMerge(const std::vector<EdgeId> &oldEdges, EdgeId newEdge) {
        int shift = 0;
        for (const auto &e : oldEdges) {
            if (edges_positions_.count(e)) {
                AddAndShiftEdgePositions(newEdge, edges_positions_[e], shift);
            }
            shift += int(this->g().length(e));
        }
    }

    virtual void HandleAdd(EdgeId /*e*/) {
    }

    virtual void HandleDelete(EdgeId e) {
        edges_positions_.erase(e);
    }

    void clear() {
        edges_positions_.clear();
    }

private:
    DECL_LOGGER("EdgesPositionHandler");
};

}

#endif /* EDGES_POSITION_HANDLER_HPP_ */
