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
            components.push_back(component);
        }
    }
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
bool ComponentEstimator::IsCorrect(const ComponentEstimator::TransitionGraph &transition_graph,
                                   const std::vector<std::set<ScaffoldVertex>> &clusters) const {
    VERIFY_DEV(transition_graph.size() >= 2);
    for (const auto &cluster: clusters) {
        if (not path_cluster_validator_.IsCorrect(cluster)) {
            return false;
        }
    }
    set<ScaffoldVertex> graph_vertices;
    std::copy(transition_graph.begin(), transition_graph.end(), std::inserter(graph_vertices, graph_vertices.end()));
    auto component_path_result = path_cluster_validator_.GetReferencePath(graph_vertices);
    if (not component_path_result.is_initialized()) {
        return false;
    }
    const auto &component_path = component_path_result.get();
    VERIFY_DEV(component_path.size() == transition_graph.size());
    for (auto first = component_path.begin(), second = std::next(first); second != component_path.end(); ++first, ++second) {
        if (not transition_graph.ContainsEdge(*first, *second)) {
            return false;
        }
    }
    return true;
}
void ComponentEstimator::EstimateComponents(const ComponentEstimator::ScaffoldGraph &scaffold_graph) const {
    ScaffoldGraphComponentExtractor component_extractor;
    auto components = component_extractor.GetConnectedComponents(scaffold_graph);
    vector<TransitionGraph> covered_components;
    vector<TransitionGraph> true_components;
    vector<TransitionGraph> simple_components;
    vector<TransitionGraph> false_components;
    vector<TransitionGraph> trivial_components;
    for (const auto &component: components) {
        if (IsTrivial(component)) {
            trivial_components.push_back(component);
            continue;
        }
        if (IsCovered(component)) {
            DEBUG("Getting final clusters");
            auto final_clusters = path_cluster_helper_.GetFinalClusters(component);
            for (const auto &cluster: final_clusters) {
                VERIFY_DEV(path_cluster_validator_.IsCovered(cluster));
            }
            if (IsCorrect(component, final_clusters)) {
                true_components.push_back(component);
                if (IsSimple(component, final_clusters)) {
                    simple_components.push_back(component);
                }
            } else {
                false_components.push_back(component);
            }
            covered_components.push_back(component);
        }
    }
    INFO("Trivial components: " << trivial_components.size());
    INFO("Covered components: " << covered_components.size());
    INFO("Correct components: " << true_components.size());
    INFO("Simple components: " << simple_components.size());
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
bool ComponentEstimator::IsSimple(const ComponentEstimator::TransitionGraph &transition_graph,
                                  const std::vector<std::set<ScaffoldVertex>> &clusters) const {
    vector<set<ScaffoldVertex>> doubletons;
    std::copy_if(clusters.begin(), clusters.end(), std::back_inserter(doubletons),
                 [this](const set<ScaffoldVertex> &cluster) {
                   return cluster.size() == 2;
                 });
    for (const auto &cluster: doubletons) {
        if (not path_cluster_validator_.IsCorrect(cluster)) {
            return false;
        }
    }
    return doubletons.size() == transition_graph.size() - 1;
}
bool ComponentEstimator::IsCovered(const ComponentEstimator::TransitionGraph &transition_graph) const {
    for (const auto &vertex: transition_graph) {
        if (not path_cluster_validator_.IsCovered(vertex)) {
            return false;
        }
    }
    return true;
}
bool ComponentEstimator::IsTrivial(const ComponentEstimator::TransitionGraph &transition_graph) const {
    return transition_graph.size() <= 2;
}
}
