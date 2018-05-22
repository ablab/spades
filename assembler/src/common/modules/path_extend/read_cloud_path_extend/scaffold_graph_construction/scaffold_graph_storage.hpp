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

        size_t large_length_threshold_;
        size_t small_length_threshold_;

     public:
        ScaffoldGraphStorage(const debruijn_graph::Graph& g);

        ScaffoldGraphStorage(ScaffoldGraph&& large_scaffold_graph, ScaffoldGraph&& small_scaffold_graph,
                             size_t large_length_threshold, size_t small_length_threshold);

        ScaffoldGraphStorage& operator= (const ScaffoldGraphStorage &other);

        const ScaffoldGraph& GetLargeScaffoldGraph() const;

        const ScaffoldGraph& GetSmallScaffoldGraph() const;

        size_t GetLargeLengthThreshold() const;

        size_t GetSmallLengthThreshold() const;

        void SetLargeScaffoldGraph(const ScaffoldGraph& large_scaffold_graph);

        void SetSmallScaffoldGraph(const ScaffoldGraph& small_scaffold_graph);

        void Save(const string& path) const;

        void Load(const string& path, const std::map<size_t, debruijn_graph::EdgeId>& edge_map);

     private:
        void ReplaceScaffoldGraph(const ScaffoldGraph &from, ScaffoldGraph &to);
    };

    class ScaffoldGraphSerializer {
     public:
        typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef debruijn_graph::EdgeId EdgeId;

        void SaveScaffoldGraph(ofstream& fout, const ScaffoldGraph& graph) const;

        void LoadScaffoldGraph(ifstream& fin, ScaffoldGraph& graph,
                               const std::map<size_t, debruijn_graph::EdgeId>& edge_map) const;

        DECL_LOGGER("ScaffoldGraphSerializer");
    };
}