//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "common/pipeline/graph_pack.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "common/modules/path_extend/path_extender.hpp"
#include <unordered_map>
#include <numeric>

namespace debruijn_graph {
    class DomainGraphConstructor;
}

namespace nrps {
    class DomainGraphDataMaster;

    class DomainVertexData {
        friend class DomainDataMaster;
        typedef debruijn_graph::EdgeId EdgeId;
        omnigraph::MappingPath<EdgeId> mapping_path_;
        std::string domain_name_;
        std::string domain_desc_;
        size_t start_coord_;
        size_t end_coord_;
        size_t max_visited_;
        size_t current_visited_;
        bool visited_;
        bool was_started_;
        bool near_contig_end_;
        bool near_contig_start_;

        omnigraph::MappingRange conjugate(omnigraph::MappingRange m, EdgeId e, const debruijn_graph::Graph &g) const {
        return omnigraph::MappingRange(m.initial_range.start_pos, m.initial_range.end_pos,
                                       g.length(e) - m.mapped_range.end_pos, g.length(e) - m.mapped_range.start_pos);
    }

    public:
        DomainVertexData()
                : domain_name_("None"),
                  start_coord_(0), end_coord_(0),
                  max_visited_(1), current_visited_(0), visited_(false), was_started_(false) {
        }

        DomainVertexData(const omnigraph::MappingPath<EdgeId> &mapping_path,
                         const std::string &domain_name,
                         const std::string &domain_desc,
                         size_t start_coord,
                         size_t end_coord,
                         size_t max_visited = 1,
                         size_t current_visited = 0,
                         bool near_contig_end = false,
                         bool near_contig_start = false,
                         bool visited = false,
                         bool was_started = false)
        : mapping_path_(mapping_path),
          domain_name_(domain_name), domain_desc_(domain_desc),
          start_coord_(start_coord), end_coord_(end_coord),
          max_visited_(max_visited), current_visited_(current_visited), visited_(visited), was_started_(was_started),
          near_contig_end_(near_contig_end), near_contig_start_(near_contig_start) {}

        size_t GetStartCoord() const { return start_coord_; }
        size_t GetEndCoord() const { return end_coord_; }

        bool GetNearContigEnd() const { return near_contig_end_; }
        bool GetNearContigStart() const { return near_contig_start_; }

        void SetNearStartCoord() { near_contig_start_ = true; }
        void SetNearEndCoord() { near_contig_end_ = true; }

        const std::string &GetDomainName() const { return domain_name_; }
        const std::string &GetDomainDesc() const { return domain_desc_; }

        size_t GetMaxVisited() const { return max_visited_; }
        void SetMaxVisited(size_t value) { max_visited_ = value; }

        void IncrementVisited() { current_visited_++; }
        void DecrementVisited() {  current_visited_--; }
        size_t GetCurrentVisited() const { return current_visited_; }

        bool Visited() const { return visited_; }
        bool WasStarted() const { return was_started_; }

        void SetVisited() { visited_ = true; }
        void SetWasStarted() { was_started_ = true; }

        const std::vector<EdgeId> &domain_edges() const {
            return mapping_path_.simple_path();
        }
        const omnigraph::MappingPath<EdgeId>& mapping_path() const {
            return mapping_path_;
        }

        DomainVertexData conjugate(const debruijn_graph::Graph &g) const {
            omnigraph::MappingPath<EdgeId> conjugate_rc;
            if (mapping_path_.size() != 0) {
                for (auto i = mapping_path_.size(); i != 0; --i) {
                    conjugate_rc.push_back(g.conjugate(mapping_path_.edge_at(i - 1)),
                                           conjugate(mapping_path_.mapping_at(i - 1),
                                                     g.conjugate(mapping_path_.edge_at(i - 1)), g));
                }
            }
            return DomainVertexData(conjugate_rc, domain_name_, domain_desc_,
                                    g.length(conjugate_rc.front().first) - end_coord_,
                                    g.length(conjugate_rc.back().first) - start_coord_, max_visited_,
                                    current_visited_, near_contig_start_, near_contig_end_, visited_, was_started_);
        }

        size_t length(const debruijn_graph::Graph &g) const {
            auto simple_path = mapping_path_.simple_path();
            return std::accumulate(simple_path.begin(), simple_path.end(), 0,
                                    [&](size_t a, EdgeId b){return a + g.length(b);}) - start_coord_ - end_coord_;
        }
    };

    class DomainEdgeData {
        friend class DomainGraphDataMaster;
        typedef debruijn_graph::EdgeId EdgeId;
        bool strong_;
        std::vector<EdgeId> edges_;
        size_t length_;
    public:

        explicit DomainEdgeData(bool strong, const std::vector<EdgeId> &edges, size_t length)
                : strong_(strong), edges_(edges), length_(length) {}

        bool strong() const { return strong_; }
        size_t length(const debruijn_graph::Graph &) const { return length_; }
        const std::vector<EdgeId> &debruijn_edges() const { return edges_;  }

        DomainEdgeData conjugate(const debruijn_graph::Graph &g) const {
            std::vector<EdgeId> rc;
            for (auto it = edges_.rbegin(); it != edges_.rend(); ++it) {
                rc.push_back(g.conjugate(*it));
            }
            return DomainEdgeData(strong_, rc, length_);
        }
    };

    class DomainGraphDataMaster {
        const debruijn_graph::Graph &g_;
    public:
        typedef DomainVertexData VertexData;
        typedef DomainEdgeData EdgeData;

        DomainGraphDataMaster(const debruijn_graph::Graph &g)
                : g_(g) { }

        EdgeData conjugate(const EdgeData &data) const {
            return data.conjugate(g_);
        }

        VertexData conjugate(const VertexData &data) const {
            return data.conjugate(g_);
        }

        size_t length(const EdgeData& data) const {
            return data.length(g_);
        }

        bool isSelfConjugate(const EdgeData &) const {
            return false;
        }

        size_t length(const VertexData& data) const {
            return data.length(g_);
        }
    };

    class DomainGraph : public omnigraph::ObservableGraph<DomainGraphDataMaster> {
      public:
        struct Arrangements {
            size_t component_size = 0;
            size_t strong_edges = 0;
            size_t weak_edges = 0;
            std::vector<std::vector<VertexId>> arrangements;

            bool empty() const {
                return component_size == 0;
            }
        };

    private:
        typedef base::VertexData VertexData;
        typedef base::EdgeData EdgeData;
        typedef omnigraph::ObservableGraph<DomainGraphDataMaster> base;

        void SetVisited(VertexId v);
        void SetWasStarted(VertexId v);
        void SetMaxVisited(VertexId v, size_t value);
        size_t GetMaxVisited(VertexId v) const;
        size_t GetCurrentVisited(VertexId v) const;
        size_t GetStartCoord(VertexId v) const;
        size_t GetEndCoord(VertexId v) const;
        const std::vector<EdgeId> &debruijn_edges(EdgeId) const;
        void IncrementVisited(VertexId v);
        void DecrementVisited(VertexId v);

        void OutputStat(const DomainGraph::Arrangements &arr, std::ostream &stat_file) const;
        size_t GetMaxVisited(VertexId v, double base_coverage) const;
        void SetCopynumber(const std::set<VertexId> &preliminary_visited);
        void OutputStatArrangement(const std::vector<VertexId> &single_candidate, unsigned id, std::ostream &stat_file);
        void FindBasicStatistic(std::ostream &stat_stream);
        void PrelimDFS(VertexId v, std::set<VertexId> &preliminary_visited);
        void PathToSequence(path_extend::BidirectionalPath &p, const std::vector<VertexId> &answer);
        DomainGraph::Arrangements FindAllPossibleArrangements(VertexId v,
                                                              size_t component_size_part, size_t component_min_size);
        void FinalDFS(VertexId v, std::vector<VertexId> &current, std::set<VertexId> preliminary_visited,
                                   std::vector<std::vector<VertexId>> &answer, size_t component_size, size_t &iteration_number);
        std::set<EdgeId> CollectEdges(const path_extend::BidirectionalPath &p) const;
        void OutputComponent(const path_extend::BidirectionalPath &p, int component_id, int ordering_id);

    public:
        DomainGraph(const debruijn_graph::Graph &g) :
                base(DomainGraphDataMaster(g)), g_(g) {}

        using base::AddVertex;
        using base::AddEdge;
        using base::EdgeId;
        using base::VertexId;

        VertexId AddVertex();
        VertexId AddVertex(const std::string &vname, const omnigraph::MappingPath<EdgeId> &mapping_path,
                           size_t start_coord, size_t end_coord,
                           const std::string &name, const std::string &desc);
        const std::string &GetVertexName(VertexId v) const;
        EdgeId AddEdge(VertexId from, VertexId to, bool strong, const std::vector<EdgeId> &edges, size_t length);
        bool HasStrongEdge(VertexId v);
        size_t StrongEdgeCount(VertexId v) const;
        size_t WeakEdgeCount(VertexId v) const;
        bool HasStrongIncomingEdge(VertexId v) const;
        bool NearContigStart(VertexId v) const;
        bool NearContigEnd(VertexId v) const;
        bool Visited(VertexId v) const;
        bool WasStarted(VertexId v) const;
        const std::string &GetDomainName(VertexId v) const;
        const std::string &GetDomainDesc(VertexId v) const;

        const std::vector<EdgeId>& domain_edges(VertexId v) const;
        const omnigraph::MappingPath<EdgeId>& mapping_path(VertexId v) const;
        bool strong(EdgeId e) const;
        void SetContigNearEnd(VertexId v);
        void ExportToDot(const std::string &output_path) const;
        void FindDomainOrderings(debruijn_graph::GraphPack &gp,
                                 size_t component_size_part, size_t component_min_size, bool start_only_from_tips,
                                 const std::string &output_filename, const std::string &output_dir);
        path_extend::Gap ConvertGapDescription(const path_extend::GapDescription &gap) const;

        friend class debruijn_graph::DomainGraphConstructor;
      protected:
        DECL_LOGGER("DomainGraph");
      private:
        const debruijn_graph::Graph &g_;
        std::unordered_map<VertexId, std::string> from_id_to_name;
        path_extend::PathContainer contig_paths_;
    };
}
