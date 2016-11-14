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

//#include "utils.hpp"
#include "utils/simple_tools.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/core/action_handlers.hpp"

namespace omnigraph {

struct EdgePosition {
    string contigId;
    MappingRange mr;
    EdgePosition(string _contigId, MappingRange _mr) : contigId(_contigId), mr(_mr) {
    }

    EdgePosition() {
    }
};

inline ostream& operator <<(ostream& os, const EdgePosition& ep) {
    return os << ep.contigId << " " << ep.mr;
}

template<class Graph>
class EdgesPositionHandler: public GraphActionHandler<Graph> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    size_t max_mapping_gap_;
    size_t max_gap_diff_;
    map<EdgeId, map<string, std::set<MappingRange>>> edges_positions_;
    //TODO extract set<MappingRange> as a storage class

    MappingRange EraseAndExtract(set<MappingRange> &ranges, set<MappingRange>::iterator &position, const MappingRange &new_pos) {
        auto &old_pos = *position;
        if(old_pos.IntersectLeftOf(new_pos) || old_pos.StrictlyContinuesWith(new_pos, max_mapping_gap_, max_gap_diff_)) {
            ranges.erase(position);
            return old_pos.Merge(new_pos);
        } else if(new_pos.IntersectLeftOf(old_pos) || new_pos.StrictlyContinuesWith(old_pos, max_mapping_gap_, max_gap_diff_)) {
            ranges.erase(position);
            return new_pos.Merge(old_pos);
        } else {
            return new_pos;
        }
    }

public:
    MappingRange EraseAndExtract(set<MappingRange> &ranges, MappingRange new_pos) {
        auto it = ranges.lower_bound(new_pos);
        if(it != ranges.end()) {
            new_pos = EraseAndExtract(ranges, it, new_pos);
            it = ranges.lower_bound(new_pos);
        }
        if(it != ranges.begin()) {
            new_pos = EraseAndExtract(ranges, --it, new_pos);
        }
        return new_pos;
    }

    set<MappingRange> GetEdgePositions(EdgeId edge, string contig_id) const {
        VERIFY(this->IsAttached());
        auto edge_it = edges_positions_.find(edge);
        if(edge_it == edges_positions_.end())
            return set<MappingRange>();
        const auto& positions = edge_it->second;
        auto it = positions.find(contig_id);
        if(it == positions.end())
            return set<MappingRange>();
        else
            return it->second;
    }

    vector<EdgePosition> GetEdgePositions(EdgeId edge) const {
        VERIFY(this->IsAttached());
        auto edge_it = edges_positions_.find(edge);
        if(edge_it == edges_positions_.end())
            return vector<EdgePosition>();
        vector<EdgePosition> result;
        for(auto it = edge_it->second.begin(); it != edge_it->second.end(); ++it) {
            for(auto pos_it = it->second.begin(); pos_it != it->second.end(); ++pos_it) {
                result.push_back(EdgePosition(it->first, *pos_it));
            }
        }
        return result;
    }

    void AddEdgePosition(EdgeId edge, string contig_id, size_t start, size_t end, size_t m_start, size_t m_end) {
        VERIFY(this->IsAttached());
        AddEdgePosition(edge, contig_id, MappingRange(start, end, m_start, m_end));
    }

    void AddEdgePosition(EdgeId edge, string contig_id, MappingRange new_pos) {
        VERIFY(this->IsAttached());
        if(new_pos.empty())
            return;
        set<MappingRange> &new_set = edges_positions_[edge][contig_id];
        new_pos = EraseAndExtract(new_set, new_pos);
        new_set.insert(new_pos);
    }

    void AddAndShiftEdgePositions(EdgeId edge, const map<string, set<MappingRange>> &contig_map, int shift = 0) {
        VERIFY(this->IsAttached());
        for(auto contig_it = contig_map.begin(); contig_it != contig_map.end(); ++contig_it) {
            for(auto it = contig_it->second.begin(); it != contig_it->second.end(); ++it) {
                AddEdgePosition(edge, contig_it->first, it->Shift(shift).Fit(this->g().length(edge)));
            }
        }
    }

    template<typename Iter>
    void AddEdgePositions(EdgeId edge, Iter begin, Iter end) {
        VERIFY(this->IsAttached());
        for(auto it = begin; it != end; ++it) {
            AddEdgePosition(edge, it->contigId, it->mr);
        }
    }

    std::string str(EdgeId edge) const {
        VERIFY(this->IsAttached());
        std::stringstream ss;
        vector<EdgePosition> positions = GetEdgePositions(edge);
        size_t counter = 0;
        for (auto pos_it = positions.begin(), end = positions.end(); pos_it != end; ++pos_it) {
            ss << "(" << pos_it->contigId << ": " << pos_it->mr << ")\\n";
            counter++;
            if(counter > 30) {
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

    virtual void HandleMerge(const vector<EdgeId>& oldEdges, EdgeId newEdge) {
        int shift = 0;
        for(auto it = oldEdges.begin(); it != oldEdges.end(); ++it) {
            if (edges_positions_.count(*it) != 0) {
                AddAndShiftEdgePositions(newEdge, edges_positions_[*it], shift);
            }
            shift += int(this->g().length(*it));
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
