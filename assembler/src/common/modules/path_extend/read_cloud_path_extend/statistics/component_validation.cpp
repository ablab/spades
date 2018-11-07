#include "component_validation.hpp"
#include <stack>
namespace path_extend {

ScaffoldGraphComponentExtractor::TransitionGraph ScaffoldGraphComponentExtractor::UnorientTransitionGraph(
        const ScaffoldGraphComponentExtractor::TransitionGraph &transition_graph) const {
    TransitionGraph result;
    for (const auto &vertex: transition_graph) {
        result.AddVertex(vertex);
    }
    for (const auto &vertex: transition_graph) {
        for (auto it = transition_graph.outcoming_begin(vertex); it != transition_graph.outcoming_end(vertex); ++it) {
            result.AddEdge(vertex, *it);
            result.AddEdge(*it, vertex);
        }
    }
    return result;
}
vector<ScaffoldGraphComponentExtractor::TransitionGraph> ScaffoldGraphComponentExtractor::GetConnectedComponents(
        const ScaffoldGraphComponentExtractor::ScaffoldGraph &scaffold_graph) const {
    TransitionGraph simple_graph;
    for (const auto &vertex: scaffold_graph.vertices()) {
        simple_graph.AddVertex(vertex);
    }
    for (const auto &edge: scaffold_graph.edges()) {
        simple_graph.AddEdge(edge.getStart(), edge.getEnd());
    }
    auto unoriented_graph = UnorientTransitionGraph(simple_graph);
    set<ScaffoldVertex> visited;
    vector<TransitionGraph> components;
    size_t trivial_components = 0;
    for (const auto &vertex: unoriented_graph) {
        if (visited.find(vertex) == visited.end()) {
            auto component_vertices = GetVertexComponent(unoriented_graph, vertex);
            std::copy(component_vertices.begin(), component_vertices.end(), std::inserter(visited, visited.end()));
            TransitionGraph component;
            DEBUG("Adding component vertices");
            for (const auto &vertex: component_vertices) {
                component.AddVertex(vertex);
            }
            DEBUG("Adding component edges");
            for (const auto &vertex: component_vertices) {
                for (auto it = simple_graph.outcoming_begin(vertex); it != simple_graph.outcoming_end(vertex); ++it) {
                    TRACE("Adding edge " << vertex.int_id() << " -> " << it->int_id());
                    component.AddEdge(vertex, *it);
                }
            }
            if (component.size() > 2) {
                components.push_back(component);
            } else {
                ++trivial_components;
            }
        }
    }
    INFO("Trivial components: " << trivial_components);
    return components;
}
set<ScaffoldGraphComponentExtractor::ScaffoldVertex> ScaffoldGraphComponentExtractor::GetVertexComponent(
    const ScaffoldGraphComponentExtractor::TransitionGraph &transition_graph,
    const ScaffoldGraphComponentExtractor::ScaffoldVertex &start) const {
    set<ScaffoldVertex> visited;
    std::stack<ScaffoldVertex> current_stack;
    current_stack.push(start);
    visited.insert(start);
    while (not current_stack.empty()) {
        ScaffoldVertex current = current_stack.top();
        current_stack.pop();
        TRACE("Processing " << current.int_id());
        for (auto it = transition_graph.outcoming_begin(current); it != transition_graph.outcoming_end(current); ++it) {
            ScaffoldVertex next = *it;
            if (visited.find(next) == visited.end()) {
                TRACE("Inserting " << next.int_id());
                visited.insert(next);
                current_stack.push(next);
            }
        }
    }
    return visited;
}
bool ComponentEstimator::IsCorrect(const ComponentEstimator::TransitionGraph &transition_graph) const {
    auto final_clusters = path_cluster_helper_.GetFinalClusters(transition_graph);
    vector<set<ScaffoldVertex>> covered_clusters;
    std::copy_if(final_clusters.begin(), final_clusters.end(), std::back_inserter(covered_clusters),
                 [this](const set<ScaffoldVertex> &cluster) {
                   return this->path_cluster_validator_.IsCovered(cluster);
                 });
    vector<set<ScaffoldVertex>> true_clusters;
    vector<set<ScaffoldVertex>> false_clusters;
    for (const auto &cluster: covered_clusters) {
        if (not path_cluster_validator_.IsCorrect(cluster)) {
            return false;
        }
    }
    return true;
}
void ComponentEstimator::EstimateComponents(const ComponentEstimator::ScaffoldGraph &scaffold_graph) const {
    ScaffoldGraphComponentExtractor component_extractor;
    auto components = component_extractor.GetConnectedComponents(scaffold_graph);
    vector<TransitionGraph> true_components;
    vector<TransitionGraph> false_components;
    for (const auto &component: components) {
        if (IsCorrect(component)) {
            true_components.push_back(component);
        } else {
            false_components.push_back(component);
        }
    }
    INFO("True components: " << true_components.size());
    INFO("False components: " << false_components.size());
    auto add_component_length = [this](size_t current_sum, const TransitionGraph &graph) {
      size_t component_length = std::accumulate(graph.begin(), graph.end(), 0,
                                                [this](size_t current_length, const ScaffoldVertex &vertex) {
                                                  return current_length + vertex.getLengthFromGraph(this->g_);
                                                });
      return current_sum + component_length;
    };
    size_t true_length = std::accumulate(true_components.begin(), true_components.end(), 0, add_component_length);
    size_t false_length = std::accumulate(false_components.begin(), false_components.end(), 0, add_component_length);
    INFO("True length: " << true_length);
    INFO("False length: " << false_length);
}
ComponentEstimator::ComponentEstimator(const Graph &g,
                                       const ScaffoldGraphPathClusterHelper &path_cluster_helper,
                                       const validation::PathClusterValidator &path_cluster_validator)
    : g_(g), path_cluster_helper_(path_cluster_helper), path_cluster_validator_(path_cluster_validator) {}
}
