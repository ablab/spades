//
// Created by lab42 on 8/24/15.
//

#include "connected_component.hpp"
#include <stack>


namespace debruijn_graph {


void ConnectedComponentCounter::CalculateComponents() const {
    map <EdgeId, size_t> component_ids;
    vector <pair<size_t, size_t>> to_sort;
    map<size_t, size_t> comp_size;
    size_t cur_id = 0;
    for (auto e = g_.ConstEdgeBegin(); !e.IsEnd(); ++e) {
        if (component_ids.find(*e) == component_ids.end()) {
            std::stack <EdgeId> next;
            next.push(*e);
            set <EdgeId> used;
            size_t ans = 0;
            while (!next.empty()) {
                auto cur = next.top();
                next.pop();
                if (used.find(cur) != used.end()) {
                    continue;
                }
                ans += g_.length(cur);
                used.insert(cur);
                vector <EdgeId> neighbours;
                neighbours.push_back(g_.conjugate(cur));
                auto start = g_.EdgeStart(cur);
                auto tmp = g_.IncidentEdges(start);

                neighbours.insert(neighbours.end(), tmp.begin(), tmp.end());
                auto end = g_.EdgeEnd(cur);
                tmp = g_.IncidentEdges(end);
                neighbours.insert(neighbours.end(), tmp.begin(), tmp.end());
                for (auto ee:neighbours) {
                    if (used.find(ee) == used.end()) {
                        next.push(ee);
                    }
                }
            }
            for (auto edge: used) {
                component_ids[edge] = cur_id;
            }
            to_sort.push_back(std::make_pair(ans, cur_id));
            comp_size[cur_id] = ans;
            cur_id ++;
        }
    }
    std::sort(to_sort.begin(), to_sort.end());
    std::reverse(to_sort.begin(), to_sort.end());
    vector <size_t> perm(to_sort.size());
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
