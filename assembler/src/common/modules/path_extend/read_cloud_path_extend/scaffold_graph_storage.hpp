#pragma once

#include "common/assembly_graph/core/graph.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include <fstream>

namespace path_extend {
    class ScaffoldGraphStorage {
     public:
        typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef debruijn_graph::EdgeId EdgeId;
     private:
        ScaffoldGraph large_scaffold_graph_;
        ScaffoldGraph small_scaffold_graph_;

     public:
        explicit ScaffoldGraphStorage(const debruijn_graph::Graph& g);
        ScaffoldGraphStorage(const ScaffoldGraph& large_scaffold_graph, const ScaffoldGraph& small_scaffold_graph);

        const ScaffoldGraph& GetLargeScaffoldGraph() const;

        const ScaffoldGraph& GetSmallScaffoldGraph() const;

        void SetLargeScaffoldGraph(const ScaffoldGraph& large_scaffold_graph);

        void SetSmallScaffoldGraph(const ScaffoldGraph& small_scaffold_graph);

        void Save(const string& path) const;

        void SaveScaffoldGraph(ofstream& fout, const ScaffoldGraph& graph) const;

        void Load(const string& path, const std::map<size_t, debruijn_graph::EdgeId>& edge_map);

        void LoadScaffoldGraph(ifstream& fin, ScaffoldGraph& graph,
                               const std::map<size_t, debruijn_graph::EdgeId>& edge_map) const;

     private:
        void CopyScaffoldGraph(const ScaffoldGraph& from, ScaffoldGraph& to);

};
}