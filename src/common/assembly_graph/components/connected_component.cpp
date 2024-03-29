
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//
// Created by lab42 on 8/24/15.
//

#include "connected_component.hpp"
#include <algorithm>
#include <map>
#include <stack>
#include <vector>
#include <set>
#include <unordered_set>

namespace debruijn_graph {

void ConnectedComponentCounter::CalculateComponents() const {
    std::map<EdgeId, size_t> component_ids;
    std::vector<std::pair<size_t, size_t>> to_sort;
    std::map<size_t, size_t> comp_size;
    size_t cur_id = 0;
    for (EdgeId e : g_.edges()) {
        if (component_ids.find(e) != component_ids.end())
            continue;
        
        std::unordered_set <EdgeId> next;
        next.insert(e);
        std::set<EdgeId> used;
        size_t ans = 0;
        while (!next.empty()) {
            auto cur = *next.begin();
            next.erase(next.begin());
            if (used.find(cur) != used.end())
                continue;

            ans += g_.length(cur);
            used.insert(cur);
            std::vector<EdgeId> neighbours;
            neighbours.push_back(g_.conjugate(cur));
            auto start = g_.EdgeStart(cur);
            auto tmp = g_.IncidentEdges(start);

            neighbours.insert(neighbours.end(), tmp.begin(), tmp.end());
            auto end = g_.EdgeEnd(cur);
            tmp = g_.IncidentEdges(end);
            neighbours.insert(neighbours.end(), tmp.begin(), tmp.end());
            for (auto ee:neighbours) {
                if (used.find(ee) == used.end()) {
                    next.insert(ee);
                }
            }
        }
        
        for (auto edge: used)
            component_ids[edge] = cur_id;

        to_sort.push_back(std::make_pair(ans, cur_id));
        comp_size[cur_id] = ans;
        cur_id ++;
    }
    
    //TODO: sort in descending order
    std::sort(to_sort.begin(), to_sort.end());
    std::reverse(to_sort.begin(), to_sort.end());
    std::vector<size_t> perm(to_sort.size());
    for (size_t i = 0; i < to_sort.size(); i++) {
        perm[to_sort[i].second] = i;
        component_total_len_[i] = comp_size[to_sort[i].second];
    }
    for (auto pair:component_ids) {
        component_ids_[pair.first] = perm[pair.second];
        component_edges_quantity_[perm[pair.second]]++;
    }
    return;
}

size_t ConnectedComponentCounter::GetComponent(EdgeId e) const {
    if (component_ids_.size() == 0) {
        CalculateComponents();
    }
    return component_ids_.at(e);
}


}
