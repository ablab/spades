//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "domain_graph.hpp"
#include "utils/filesystem/path_helper.hpp"
namespace nrps {
    void DomainGraph::SetVisited(VertexId v) {
        data(v).SetVisited();
        data(conjugate(v)).SetVisited();
    }

    void DomainGraph::SetWasStarted(VertexId v) {
        data(v).SetWasStarted();
        data(conjugate(v)).SetWasStarted();
    }

    void DomainGraph::SetMaxVisited(VertexId v, size_t value) {
        data(v).SetMaxVisited(value);
    }

    size_t DomainGraph::GetMaxVisited(VertexId v) const {
        return data(v).GetMaxVisited();
    }

    const std::string &DomainGraph::GetDomainName(VertexId v) const {
        return data(v).GetDomainName();
    }
    const std::string &DomainGraph::GetDomainDesc(VertexId v) const {
        return data(v).GetDomainDesc();
    }

    bool DomainGraph::HasStrongEdge(VertexId v) {
        for (auto e : OutgoingEdges(v)) {
            if (strong(e))
                return true;
        }

        return false;
    }

    size_t DomainGraph::StrongEdgeCount(VertexId v) const {
        return std::count_if(this->OutgoingEdges(v).begin(), OutgoingEdges(v).end(),
                             [&](EdgeId e) { return strong(e); });
    }

    size_t DomainGraph::WeakEdgeCount(VertexId v) const {
        return std::count_if(this->OutgoingEdges(v).begin(), OutgoingEdges(v).end(),
                             [&](EdgeId e) { return !strong(e); });
    }

    bool DomainGraph::HasStrongIncomingEdge(VertexId v) const {
        for (auto e : IncomingEdges(v)) {
            if (strong(e))
                return true;
        }

        return false;
    }

    bool DomainGraph::NearContigStart(VertexId v) const {
        return data(v).GetNearContigStart();
    }

    bool DomainGraph::NearContigEnd(VertexId v) const {
        return data(v).GetNearContigEnd();
    }

    bool DomainGraph::Visited(VertexId v) const {
        return data(v).Visited();
    }

    bool DomainGraph::WasStarted(VertexId v) const {
        return data(v).WasStarted();
    }

    const std::vector<DomainGraph::EdgeId> &DomainGraph::domain_edges(VertexId v) const {
        return data(v).domain_edges();
    }

    const omnigraph::MappingPath<debruijn_graph::EdgeId>& DomainGraph::mapping_path(VertexId v) const {
        return data(v).mapping_path();
    }

    bool DomainGraph::strong(EdgeId e) const {
        return data(e).strong();
    }

    void DomainGraph::SetContigNearEnd(VertexId v) {
        data(v).SetNearEndCoord();
        data(conjugate(v)).SetNearStartCoord();
    }

    size_t DomainGraph::GetCurrentVisited(VertexId v) const {
        return data(v).GetCurrentVisited();
    }

    size_t DomainGraph::GetStartCoord(VertexId v) const {
        return data(v).GetStartCoord();
    }

    size_t DomainGraph::GetEndCoord(VertexId v) const {
        return data(v).GetEndCoord();
    }

    const std::vector<debruijn_graph::EdgeId> &DomainGraph::debruijn_edges(EdgeId e) const {
        return data(e).debruijn_edges();
    }

    void DomainGraph::IncrementVisited(VertexId v) {
        data(v).IncrementVisited();
    }

    void DomainGraph::DecrementVisited(VertexId v) {
        data(v).DecrementVisited();
    }

    size_t DomainGraph::GetMaxVisited(VertexId v, double base_coverage) const  {
        double low_coverage = std::numeric_limits<double>::max();
        for (EdgeId e : domain_edges(v))
            low_coverage = std::min(low_coverage, g_.coverage(e));

        return size_t(round(low_coverage / base_coverage));
    }

    void DomainGraph::SetCopynumber(const std::set<VertexId> &preliminary_visited) {
        double base_coverage = std::numeric_limits<double>::max();
        for (VertexId v : preliminary_visited) {
            for (EdgeId e : this->domain_edges(v)) {
                if (g_.length(e) <= 500)
                    continue;

                base_coverage = std::min(base_coverage, g_.coverage(e));
            }
        }

        if (math::eq(base_coverage, std::numeric_limits<double>::max()))
            return;

        for (VertexId v : preliminary_visited) {
            size_t value = std::max(size_t(1), GetMaxVisited(v, base_coverage));
            if (GetMaxVisited(v) != value) {
                DEBUG(GetVertexName(v) << " copynumber has changed from " << GetMaxVisited(v)
                      << " to " << value);
            }
            SetMaxVisited(v, value);
        }
    }

    void DomainGraph::OutputStatArrangement(const std::vector<VertexId> &single_candidate,
                                            unsigned id, std::ostream &stat_file) {
        stat_file << "BGC candidate " << id << std::endl;
        std::string delimeter = "";
        bool is_nrps = false;
        bool is_pks = false;
        for (VertexId v : single_candidate) {
            stat_file << delimeter;
            delimeter = "-";
            const std::string &name = GetDomainName(v);
            if (name == "AMP") {
                stat_file << "A";
                is_nrps = true;
            } else if (name == "CStart") {
                stat_file << "C";
                is_nrps = true;
            } else if (name == "AT") {
                stat_file << "AT";
                is_pks = true;
            } else if (name == "TE") {
                stat_file << "TE";
            } else if (name == "KR") {
                stat_file << "KR";
                is_pks = true;
            } else if (name == "KS") {
                stat_file << "KS";
                is_pks = true;
            } else {
                stat_file << name;
            }
        }
        stat_file << std::endl;
        stat_file << "Predicted type: ";
        if (is_nrps && is_pks) {
            stat_file << "NRPS/PKS";
        } else if (is_nrps) {
            stat_file << "NRPS";
        } else if (is_pks) {
            stat_file << "PKS";
        } else {
            stat_file << "Custom";
        }
        stat_file << std::endl;
        stat_file << "Domain coordinates:" << std::endl;
        size_t current_coord = 0;
        auto p = path_extend::BidirectionalPath::create(g_);
        for (size_t i = 0; i < single_candidate.size(); ++i) {
            auto v = single_candidate[i];
            DEBUG("Translating vertex " << GetVertexName(v));
            for (EdgeId e : domain_edges(v)) {
                if (p->Size() == 0 || p->Back() != e) {
                    int gap = 0;
                    if (p->Size() != 0 &&
                        p->graph().EdgeEnd(p->Back()) != g_.EdgeStart(e)) {
                        gap = 100;
                    }
                    p->PushBack(e, path_extend::Gap(gap));
                    current_coord += g_.length(e) + gap;
                }
            }

            size_t sum = 0;
            for (EdgeId e2 : domain_edges(v)) {
                sum += g_.length(e2);
            }

            stat_file << GetStartCoord(v) + current_coord - sum << " ";
            stat_file << GetEndCoord(v) + current_coord - g_.length(p->Back()) << std::endl;

            //TODO: check if can be improved
            if (i != single_candidate.size() - 1) {
                auto next_edges = this->GetEdgesBetween(v, single_candidate[i + 1]);
                if (next_edges.size() == 0) {
                    continue;
                }
                EdgeId next_edge = next_edges[0];
                for (auto e : debruijn_edges(next_edge)) {
                    if (p->Size() == 0 || p->Back() != e) {
                        int gap = 0;
                        if (p->Size() != 0 &&
                            g_.EdgeEnd(p->Back()) != g_.EdgeStart(e)) {
                            gap = 100;
                        }
                        p->PushBack(e, path_extend::Gap(gap));
                        current_coord += g_.length(e) + gap;
                    }
                }
            }
        }
    }

    void DomainGraph::FindBasicStatistic(std::ostream &stat_stream) {
        // FIXME: ugly, have common source of domain information!
        std::map<std::string, std::string> domains;
        for (VertexId v : vertices())
            domains[GetDomainName(v)] = GetDomainDesc(v);

        for (const auto &entry : domains)
            stat_stream << entry.first << " - " << entry.second << std::endl;

        std::map<std::string, unsigned> domain_count;
        for (VertexId v : vertices()) {
            domain_count[GetDomainName(v)] += 1;
        }

        for (auto domain_type : domain_count) {
            stat_stream << "# " << domain_type.first << "-domains - "
                        << domain_type.second << std::endl;
        }
        stat_stream << std::endl;
    }

    void DomainGraph::PrelimDFS(VertexId v, std::set<VertexId> &preliminary_visited) {
        for (EdgeId e : OutgoingEdges(v)) {
            VertexId to = EdgeEnd(e);
            auto inserted = preliminary_visited.insert(to);
            if (inserted.second)
                PrelimDFS(to, preliminary_visited);
        }
    }

    path_extend::Gap DomainGraph::ConvertGapDescription(const path_extend::GapDescription &gap) const {
        if (gap == path_extend::GapDescription()) {
            return path_extend::Gap::INVALID();
        }
        return path_extend::Gap(gap.estimated_dist() + int(g_.k())
                   - int(gap.left_trim()) - int(gap.right_trim()),
                   {uint32_t(gap.left_trim()), uint32_t(gap.right_trim())}, false);
    }

    void DomainGraph::PathToSequence(path_extend::BidirectionalPath &p,
                                     const std::vector<VertexId> &answer) {
        auto &scaf_params = cfg::get().pe_params.param_set.scaffolder_options;
        path_extend::LAGapAnalyzer gap_analyzer(p.g(), scaf_params.min_overlap_length, scaf_params.flank_multiplication_coefficient,
                                                scaf_params.flank_addition_coefficient);
        DEBUG("Translating " << p.GetId() << " to sequence");
        for (size_t i = 0; i < answer.size(); ++i) {
            auto v = answer[i];
            DEBUG("Translating vertex " << GetVertexName(v));
            DEBUG(domain_edges(v));
            for (EdgeId e : domain_edges(v)) {
                if (p.Size() == 0 || p.Back() != e) {
                    path_extend::Gap gap;
                    if (p.Size() != 0 && g_.IsDeadEnd(g_.EdgeEnd(p.Back())) && g_.IsDeadStart(g_.EdgeStart(e)) &&
                        g_.EdgeEnd(p.Back()) != g_.EdgeStart(e)) {
                        DEBUG("Fixing gap between " << p.Back() << " and " << e);
                        omnigraph::GapDescription<debruijn_graph::Graph> gap_description(p.Back(), e, -(int)g_.k());
                        gap_analyzer.FixGap(gap_description);
                        gap = ConvertGapDescription(gap_description);
                    }
                    p.PushBack(e, gap);
                }
            }

            if (i != answer.size() - 1) {
                std::vector<EdgeId> next_edges = GetEdgesBetween(v, answer[i + 1]);
                DEBUG("Translating edge between " << GetVertexName(v) << " and " << GetVertexName(answer[i + 1]));
                if (next_edges.size() == 0) {
                    DEBUG("Something strange!");
                    continue;
                }

                EdgeId next_edge = next_edges[0];
                DEBUG(debruijn_edges(next_edge));
                for (debruijn_graph::EdgeId de : debruijn_edges(next_edge)) {
                    if (p.Size() == 0 || p.Back() != de) {
                        path_extend::Gap gap;
                        if (p.Size() != 0 &&
                            g_.EdgeEnd(p.Back()) != g_.EdgeStart(de)) {
                            DEBUG("Fixing gap between " << p.Back() << " and " << de);
                            omnigraph::GapDescription<debruijn_graph::Graph> gap_description(p.Back(), de, -(int)g_.k());
                            gap_analyzer.FixGap(gap_description);
                            gap = ConvertGapDescription(gap_description);
                        }
                        p.PushBack(de, gap);
                    }
                }
            }
        }
    }


DomainGraph::Arrangements DomainGraph::FindAllPossibleArrangements(VertexId v,
                                                                   size_t component_size_part, size_t component_min_size) {
        DEBUG("Starting from " << GetVertexName(v));
        std::set<VertexId> preliminary_visited;
        preliminary_visited.insert(v);
        PrelimDFS(v, preliminary_visited);
        size_t component_size = preliminary_visited.size();
        DEBUG("Component size " << component_size);
        if (component_size < component_min_size)
            return {};

        DomainGraph::Arrangements result;
        result.component_size = component_size;
        for (VertexId v : preliminary_visited) {
            result.strong_edges += StrongEdgeCount(v);
            result.weak_edges += WeakEdgeCount(v);
        }

        if (cfg::get().hm->set_copynumber)
            SetCopynumber(preliminary_visited);

        for (VertexId v : preliminary_visited)
            SetVisited(v);

        preliminary_visited.clear();
        std::vector<VertexId> current;
        size_t iteration_number = 0;
        FinalDFS(v, current, preliminary_visited, result.arrangements, component_size / component_size_part, iteration_number);

        return result;
    }

    void DomainGraph::FinalDFS(VertexId v, std::vector<VertexId> &current,
                               std::set<VertexId> preliminary_visited,
                               std::vector<std::vector<VertexId>> &answer,
                               size_t accepted_component_size, size_t &iteration_number) {
        iteration_number++;
        current.push_back(v);
        IncrementVisited(v);
        if (answer.size() > 50 || iteration_number > 10000000)
            return;

        bool was_extended = false;
        if (HasStrongEdge(v)) {
            for (EdgeId e : OutgoingEdges(v)) {
                if (strong(e)) {
                    VertexId to = EdgeEnd(e);
                    if (preliminary_visited.count(to))
                        continue;

                    was_extended = true;
                    FinalDFS(to, current, preliminary_visited, answer, accepted_component_size, iteration_number);
                }
            }
        } else {
            for (EdgeId e : OutgoingEdges(v)) {
                VertexId to = EdgeEnd(e);
                if (GetCurrentVisited(to) >= GetMaxVisited(to))
                    continue;

                if (preliminary_visited.count(to))
                    continue;

                was_extended = true;
                FinalDFS(to, current, preliminary_visited, answer, accepted_component_size, iteration_number);
            }
        }

        if (current.size() >= accepted_component_size && !was_extended && (current.size() == 1 || !WasStarted(current.back())))
            answer.push_back(current);

        DecrementVisited(v);
        current.pop_back();
    }

    std::set<DomainGraph::EdgeId> DomainGraph::CollectEdges(const path_extend::BidirectionalPath &p) const {
        std::set<EdgeId> edge_set;
        for (size_t i = 0; i < p.Size(); ++i) {
            edge_set.insert(p.At(i));
        }
        return edge_set;
    }

    void DomainGraph::OutputComponent(const path_extend::BidirectionalPath &p, int component_id,
                                      int ordering_id) {
        auto edges = CollectEdges(p);
        auto comp = omnigraph::GraphComponent<debruijn_graph::Graph>::FromEdges(g_,
                                                                                edges.begin(), edges.end(), true);
        std::ofstream os(cfg::get().output_dir + "/bgc_in_gfa/" +
                         std::to_string(component_id) + "_" + std::to_string(ordering_id) + ".gfa");
        gfa::GFAComponentWriter writer(comp, os);
        writer.WriteSegmentsAndLinks();
    }

    DomainGraph::VertexId DomainGraph::AddVertex() {
        return AddVertex(VertexData());
    }

    DomainGraph::VertexId DomainGraph::AddVertex(const std::string &vname, const omnigraph::MappingPath<EdgeId> &mapping_path,
                                                 size_t start_coord, size_t end_coord,
                                                 const std::string &name, const std::string &desc) {
        auto v = AddVertex(VertexData(mapping_path, name, desc, start_coord, end_coord));
        from_id_to_name[v] = vname;
        from_id_to_name[conjugate(v)] = vname + "_rc";
        return v;
    }

    const std::string &DomainGraph::GetVertexName(DomainGraph::VertexId v) const {
        return this->from_id_to_name.at(v);
    }

    DomainGraph::EdgeId DomainGraph::AddEdge(VertexId from, VertexId to, bool strong, const std::vector<EdgeId> &edges, size_t length) {
        return AddEdge(from, to, EdgeData(strong, edges, length));
    }

    void DomainGraph::ExportToDot(const std::string &output_path) const {
        std::ofstream out(output_path);
        out << "digraph domain_graph {" << std::endl;
        for (VertexId v : vertices()) {
            out << "\"" << GetVertexName(v) << "\""
                << " [label=\"" << GetDomainName(v) << " " << GetVertexName(v) << " "
                << GetMaxVisited(v) << "\"];" << std::endl;
        }

        for (EdgeId e : edges()) {
            out << "\"" << GetVertexName(EdgeStart(e)) << "\""
                << " -> "
                << "\"" << GetVertexName(EdgeEnd(e)) << "\""
                << " [label="
                << length(e)
                << " style=" << (strong(e) ? "bold" : "dotted")
                << "];" << std::endl;
        }
        out << "}";
    }


void DomainGraph::OutputStat(const DomainGraph::Arrangements &arr, std::ostream &stat_file) const {
    stat_file << "# domains in the component - "
              << arr.component_size << std::endl
              << "# Strong/weak edges in the component - "
              << arr.strong_edges << "/" << arr.weak_edges << std::endl;
}

void DomainGraph::FindDomainOrderings(debruijn_graph::GraphPack &gp,
                                      size_t component_size_part, size_t component_min_size, bool start_only_from_tips,
                                      const std::string &output_filename, const std::string &output_dir) {
    const auto &graph = gp.get<debruijn_graph::Graph>();
    std::ofstream stat_stream(fs::append_path(output_dir, "bgc_statistics.txt"));
    FindBasicStatistic(stat_stream);
    path_extend::ContigWriter writer(graph, path_extend::MakeContigNameGenerator(cfg::get().mode, gp));

    std::vector<VertexId> start_nodes;
    std::vector<VertexId> additional_start_nodes;

    for (VertexId v : vertices()) {
        if (IsDeadStart(v) &&
            !(component_min_size > 1 && !IsDeadEnd(v)))
            start_nodes.push_back(v);
    }
    if (!start_only_from_tips) {
        for (VertexId v : vertices()) {
            if ((!HasStrongIncomingEdge(v) && IncomingEdgeCount(v)))
                additional_start_nodes.push_back(v);
        }
    }

    INFO("Total start vertices: " << start_nodes.size());

    io::osequencestream_bgc oss(fs::append_path(output_dir, output_filename));
    path_extend::ScaffoldSequenceMaker seq_maker(graph);
    std::vector<DomainGraph::Arrangements> answer;

    for (VertexId v : start_nodes) {
        if (WasStarted(v))
            continue;

        auto res = FindAllPossibleArrangements(v,
                                               component_size_part, component_min_size);
        SetWasStarted(v);
        if (!res.empty())
            answer.emplace_back(std::move(res));
    }

    for (VertexId v : additional_start_nodes) {
        if (Visited(v))
            continue;

        auto res = FindAllPossibleArrangements(v,
                                               component_size_part, component_min_size);
        if (!res.empty())
            answer.emplace_back(std::move(res));
    }

    std::stable_sort(answer.begin(), answer.end(),
                     [](const DomainGraph::Arrangements &lhs, const DomainGraph::Arrangements &rhs) {
                         return lhs.component_size > rhs.component_size;
                     });

    unsigned ordering_id = 1;
    unsigned component_id = 1;

    path_extend::PathWriter path_writer(graph);
    for (const auto &entry : answer) {
        stat_stream << "BGC subgraph " << component_id << std::endl;
        OutputStat(entry, stat_stream);

        ordering_id = 1;
        for (const auto &vec : entry.arrangements) {
            OutputStatArrangement(vec, ordering_id, stat_stream);
            auto &p = contig_paths_.Create(graph);
            PathToSequence(p, vec);
            stat_stream << "Edge order: \n" << path_writer.ToPathString(p) << "\n"
                        << "Path is " << (p.IsCircular() ? "circular" : "linear") << "\n";
            OutputComponent(p, component_id, ordering_id);
            std::string outputstring = seq_maker.MakeSequence(p);
            oss.SetCluster(component_id, ordering_id, unsigned(vec.size()));
            oss << outputstring;
            ordering_id++;
        }
        component_id++;
    }
}

}
