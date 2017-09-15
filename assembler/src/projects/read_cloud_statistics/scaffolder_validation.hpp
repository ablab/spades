#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "scaffold_graph_utils.hpp"

namespace scaffold_graph_utils {

struct NamedScaffoldGraph {
  string name_;
  ScaffoldGraph graph;
  NamedScaffoldGraph(const string& name_, const ScaffoldGraph& graph) : name_(name_), graph(graph) {}
};

class ScaffoldGraphStorage {
 public:
    typedef std::map<string, NamedScaffoldGraph>::const_iterator const_iterator;
 protected:
    std::map<string, NamedScaffoldGraph> name_to_graph_;

 public:
    void insert(const string& name, const ScaffoldGraph& graph) {
        NamedScaffoldGraph named_graph(name, graph);
        name_to_graph_.insert({name, named_graph});
    }

    NamedScaffoldGraph at(const string& name) const {
        return name_to_graph_.at(name);
    }

    const_iterator begin() const {
        return name_to_graph_.begin();
    }

    const_iterator end() const {
        return name_to_graph_.end();
    }

    adt::iterator_range<const_iterator> entries() const {
        return adt::make_range(begin(), end());
    };
};


class ScaffolderStats: public read_cloud_statistics::Statistic {
    std::map<string, path_extend::validation::ScaffoldGraphStats> name_to_stats_;
    size_t reference_transitions_;
 public:
    ScaffolderStats(const map<string, path_extend::validation::ScaffoldGraphStats>& name_to_stats_, size_t reference_transitions)
        : Statistic("scaffolder_stats"), name_to_stats_(name_to_stats_),
          reference_transitions_(reference_transitions) {}

    void Serialize(const string& path) override {
        ofstream fout(path);
        fout << "Reference transitions: " << reference_transitions_ << endl << endl;
        for (const auto& entry: name_to_stats_) {
            fout << entry.first << std::endl << std::endl;
            entry.second.Serialize(fout);
            fout << endl << endl;
        }
    }
};

class ScaffolderAnalyzer : public read_cloud_statistics::StatisticProcessor {
 public:
    typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
    typedef path_extend::validation::ContigTransitionStorage ContigTransitionStorage;
 private:
    vector<vector<EdgeWithMapping>> reference_paths_;
    ScaffoldGraphStorage storage_;
    const Graph& g_;
 public:
    ScaffolderAnalyzer(const vector<vector<EdgeWithMapping>>& reference_paths_,
                       const ScaffoldGraphStorage& storage, const Graph& g)
        : StatisticProcessor("scaffolder_analyzer"), reference_paths_(reference_paths_), storage_(storage), g_(g) {}

    void FillStatistics() override {
        auto scaffolder_stats = std::make_shared<ScaffolderStats>(GetScaffolderStats(storage_));
        AddStatistic(scaffolder_stats);
    }

 private:

    ScaffolderStats GetScaffolderStats(const ScaffoldGraphStorage& storage) {
        path_extend::validation::ScaffoldGraphValidator validator(g_);
        INFO("Getting name to stats");
        std::map<string, path_extend::validation::ScaffoldGraphStats> name_to_stats;
        for (const auto& entry: storage.entries()) {
            string name = entry.first;
            INFO("Getting stats for " << name);
            name_to_stats.insert({name, validator.GetScaffoldGraphStats(entry.second.graph, reference_paths_)});
        }
        path_extend::validation::GeneralTransitionStorageBuilder reference_transition_builder(g_, 1, false, false);
        auto reference_transitions = reference_transition_builder.GetTransitionStorage(reference_paths_);
        ScaffolderStats result(name_to_stats, reference_transitions.size());
        for (const auto& path: reference_paths_) {
            for (const auto& edge: path) {
                std::cout << edge.edge_.int_id() << std::endl;
            }
        }
        return result;
    };

    DECL_LOGGER("ScaffolderAnalyzer");
};
}