//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "domain_graph.hpp"
namespace nrps {
    void DomainGraph::SetVisited(VertexId v) {
        this->data(v).SetVisited();
        this->data(this->conjugate(v)).SetVisited();
    }

    void DomainGraph::OutputStat(std::set<VertexId> &preliminary_visited, std::ofstream &stat_file) const {
        stat_file << "# domains in the component - ";
        stat_file << preliminary_visited.size() << std::endl;
        int strong_edge_count =
                std::accumulate(preliminary_visited.begin(), preliminary_visited.end(),
                                0, [&](int prev, VertexId v) {
                            return prev + (this->StrongEdgeCount(v));
                        });
        int weak_edge_count =
                std::accumulate(preliminary_visited.begin(), preliminary_visited.end(),
                                0, [&](int prev, VertexId v) {
                            return prev + (this->WeakEdgeCount(v));
                        });

        stat_file << "# Strong/weak edges in the component - ";
        stat_file << strong_edge_count << "/" << weak_edge_count << std::endl;
    }

    size_t DomainGraph::GetMaxVisited(VertexId v, double base_coverage) const  {
        double low_coverage = std::numeric_limits<double>::max();
        for (auto e : this->GetDomainEdges(v))
            low_coverage = std::min(low_coverage, g_.coverage(e));

        return size_t(round(low_coverage / base_coverage));
    }

    void DomainGraph::SetCopynumber(std::set<VertexId> &preliminary_visited) {
        double base_coverage = std::numeric_limits<double>::max();
        for (auto v : preliminary_visited) {
            for (auto e : this->GetDomainEdges(v)) {
                if (g_.length(e) > 500) {
                    base_coverage = std::min(base_coverage, g_.coverage(e));
                }
            }
        }

        if (math::eq(base_coverage, std::numeric_limits<double>::max()))
            return;

        for (auto v : preliminary_visited) {
            size_t value = std::max(size_t(1), GetMaxVisited(v, base_coverage));
            if (this->data(v).GetMaxVisited() != value) {
                DEBUG(this->GetVertexName(v) << " copynumber has changed from " << this->data(v).GetMaxVisited()
                                             << " to " << value);
            }
            this->data(v).SetMaxVisited(value);
        }
    }

    void DomainGraph::OutputStatArrangement(std::vector<VertexId> single_candidate,
                               int id, std::ofstream &stat_file) {
        stat_file << "BGC candidate " << id << std::endl;
        std::string delimeter = "";
        bool is_nrps = false;
        bool is_pks = false;
        for (auto v : single_candidate) {
            stat_file << delimeter;
            delimeter = "-";
            std::string type = this->data(v).GetType();
            if (type == "AMP") {
                stat_file << "A";
                is_nrps = true;
            } else if (type == "CStart") {
                stat_file << "C";
                is_nrps = true;
            } else if (type == "AT") {
                stat_file << "AT";
                is_pks = true;
            } else if (type == "TE") {
                stat_file << "TE";
            } else if (type == "KR") {
                stat_file << "KR";
                is_pks = true;
            } else if (type == "KS") {
                stat_file << "KS";
                is_pks = true;
            } else {
                stat_file << type;
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
        stat_file << "Domain cordinates:" << std::endl;
        size_t current_coord = 0;
        path_extend::BidirectionalPath p(g_);
        for (size_t i = 0; i < single_candidate.size(); ++i) {
            auto v = single_candidate[i];
            DEBUG("Translating vertex " << this->GetVertexName(v));
            for (auto e : this->GetDomainEdges(v)) {
                if (p.Size() == 0 || p.Back() != e) {
                    int gap = 0;
                    if (p.Size() != 0 &&
                        p.graph().EdgeEnd(p.Back()) != g_.EdgeStart(e)) {
                        gap = 100;
                    }
                    p.PushBack(e, path_extend::Gap(gap));
                    current_coord += g_.length(e) + gap;
                }
            }

            size_t sum = 0;
            for (auto e2 : this->GetDomainEdges(v)) {
                sum += g_.length(e2);
            }

            stat_file << this->data(v).GetStartCoord() + current_coord - sum << " ";
            stat_file << this->data(v).GetEndCoord() + current_coord - g_.length(p.Back()) << std::endl;
            //TODO: check if can be improved
            if (i != single_candidate.size() - 1) {
                auto next_edges = this->GetEdgesBetween(v, single_candidate[i + 1]);
                if (next_edges.size() == 0) {
                    continue;
                }
                EdgeId next_edge = next_edges[0];
                for (auto e : this->data(next_edge).GetDeBruijnEdges()) {
                    if (p.Size() == 0 || p.Back() != e) {
                        int gap = 0;
                        if (p.Size() != 0 &&
                            g_.EdgeEnd(p.Back()) != g_.EdgeStart(e)) {
                            gap = 100;
                        }
                        p.PushBack(e, path_extend::Gap(gap));
                        current_coord += g_.length(e) + gap;
                    }
                }
            }
        }
    }

    void DomainGraph::FindBasicStatistic(std::ofstream &stat_stream) {
        stat_stream << "A - Adenylation domain" << std::endl;
        stat_stream << "AT - Acyltransferase domain" << std::endl;
        stat_stream << "C - Condensation domain" << std::endl;
        stat_stream << "KR - Keto-reductase domain" << std::endl;
        stat_stream << "KS - Keto-synthase domain" << std::endl;
        stat_stream << "TE - Termination domain" << std::endl;
        stat_stream << std::endl;

        std::map<std::string, int> domain_count;
        for (auto iter = this->SmartVertexBegin(); ! iter.IsEnd(); ++iter)
            domain_count[this->data(*iter).GetType()]++;

        for (auto domain_type : domain_count) {
            stat_stream << "# " << domain_type.first << "-domains - "
                        << domain_type.second << std::endl;
        }
        stat_stream << std::endl;
    }

    void DomainGraph::PrelimDFS(VertexId v, std::set<VertexId> &preliminary_visited) {
        for (auto e : this->OutgoingEdges(v)) {
            VertexId to = this->EdgeEnd(e);
            if (preliminary_visited.find(to) == preliminary_visited.end()) {
                preliminary_visited.insert(to);
                PrelimDFS(to, preliminary_visited);
            }
        }
    }

    std::string DomainGraph::PathToSequence(path_extend::BidirectionalPath *p,
                               std::vector<VertexId> &answer) {
        std::stringstream ss;
        DEBUG("Translating " << p->GetId() << " to sequnece");
        for (size_t i = 0; i < answer.size(); ++i) {
            auto v = answer[i];
            DEBUG("Translating vertex " << this->GetVertexName(v));
            p->PrintDEBUG();
            DEBUG(this->GetDomainEdges(v));
            for (auto e : this->GetDomainEdges(v)) {
                if (p->Size() == 0 || p->Back() != e) {
                    int gap = 0;
                    if (p->Size() != 0 &&
                        g_.EdgeEnd(p->Back()) != g_.EdgeStart(e)) {
                        gap = 100;
                    }
                    p->PushBack(e, path_extend::Gap(gap));
                }
            }
            if (i != answer.size() - 1) {
                std::vector<EdgeId> next_edges = this->GetEdgesBetween(v, answer[i + 1]);
                DEBUG("Translating edge between " <<
                                                  this->GetVertexName(v) << " and " << this->GetVertexName(answer[i + 1]));
                if (next_edges.size() == 0) {
                    DEBUG("Something strange!");
                    continue;
                }

                EdgeId next_edge = next_edges[0];
                for (auto e : this->data(next_edge).GetDeBruijnEdges()) {
                    if (p->Size() == 0 || p->Back() != e) {
                        int gap = 0;
                        if (p->Size() != 0 &&
                            g_.EdgeEnd(p->Back()) != g_.EdgeStart(e)) {
                            gap = 100;
                        }
                        p->PushBack(e, path_extend::Gap(gap));
                    }
                }
            }
        }
        return ss.str();
    }

    void DomainGraph::FindAllPossibleArrangements(VertexId v, std::vector<std::vector<VertexId>> &answer,
                                     std::ofstream &stat_file) {
        DEBUG("Starting from " << this->GetVertexName(v));
        std::set<VertexId> preliminary_visited;
        preliminary_visited.insert(v);
        PrelimDFS(v, preliminary_visited);
        size_t component_size = preliminary_visited.size();
        DEBUG("Component size " << component_size);
        if (component_size == 1)
            return;

        SetCopynumber(preliminary_visited);
        OutputStat(preliminary_visited, stat_file);
        std::vector<VertexId> current;
        for (auto v : preliminary_visited) {
            SetVisited(v);
        }
        preliminary_visited.clear();
        size_t iteration_number = 0;
        FinalDFS(v, current, preliminary_visited, answer, component_size, iteration_number);
    }

    void DomainGraph::FinalDFS(VertexId v, std::vector<VertexId> &current,
                  std::set<VertexId> preliminary_visited,
                  std::vector<std::vector<VertexId>> &answer,
                  size_t component_size, size_t &iteration_number) {
        iteration_number++;
        current.push_back(v);
        this->data(v).IncrementVisited();
        if (answer.size() > 50 || iteration_number > 10000000)
            return;

        bool was_extended = false;
        if (this->HasStrongEdge(v)) {
            for (EdgeId e : this->OutgoingEdges(v)) {
                if (this->Strong(e)) {
                    VertexId to = this->EdgeEnd(e);
                    if (preliminary_visited.find(to) != preliminary_visited.end())
                        continue;
                    was_extended = true;
                    FinalDFS(to, current, preliminary_visited, answer, component_size, iteration_number);
                }
            }
        } else {
            for (EdgeId e : this->OutgoingEdges(v)) {
                VertexId to = this->EdgeEnd(e);
                if (this->data(to).GetCurrentVisited() >= this->data(to).GetMaxVisited())
                    continue;

                if (preliminary_visited.find(to) != preliminary_visited.end())
                    continue;

                was_extended = true;
                FinalDFS(to, current, preliminary_visited, answer, component_size, iteration_number);
            }
        }
        if (current.size() >= component_size && was_extended == false) {
            answer.push_back(current);
        }
        this->data(v).DecrementVisited();
        current.pop_back();
    }

    std::set<DomainGraph::EdgeId> DomainGraph::CollectEdges(path_extend::BidirectionalPath *p) const {
        std::set<EdgeId> edge_set;
        for (size_t i = 0; i < p->Size(); ++i) {
            edge_set.insert(p->At(i));
        }
        return edge_set;
    }

    void DomainGraph::OutputComponent(path_extend::BidirectionalPath *p, int component_id,
                         int ordering_id) {
        auto edges = CollectEdges(p);
        auto comp = omnigraph::GraphComponent<debruijn_graph::Graph>::FromEdges(g_, edges.begin(),
                                                                     edges.end(), true);
        std::ofstream os(cfg::get().output_dir + "/bgc_in_gfa/" +
                         std::to_string(component_id) + "_" + std::to_string(ordering_id) +
                         ".gfa");
        gfa::GFAComponentWriter writer(comp, os);
        writer.WriteSegmentsAndLinks();
    }

    DomainGraph::VertexId DomainGraph::AddVertex() {
        return AddVertex(VertexData(g_));
    }

    DomainGraph::VertexId DomainGraph::AddVertex(const std::string &name, const omnigraph::MappingPath<EdgeId> &mapping_path,
                       size_t start_coord, size_t end_coord, std::string type) {
        auto v = AddVertex(VertexData(g_, mapping_path, type, start_coord, end_coord));
        from_id_to_name[v] = name;
        from_id_to_name[conjugate(v)] = name + "_rc";
        return v;
    }

    std::string DomainGraph::GetVertexName(DomainGraph::VertexId v) const {
        return this->from_id_to_name.at(v);
    }

    DomainGraph::EdgeId DomainGraph::AddEdge(VertexId from, VertexId to, bool strong, const std::vector<EdgeId> &edges, size_t length) {
        return AddEdge(from, to, EdgeData(g_, strong, edges, length));
    }

    bool DomainGraph::HasStrongEdge(VertexId v) {
        auto edges = this->OutgoingEdges(v);
        for (auto e : edges) {
            if (this->Strong(e)) {
                return true;
            }
        }
        return false;
    }

    size_t DomainGraph::StrongEdgeCount(VertexId v) const {
        return std::count_if(this->OutgoingEdges(v).begin(), this->OutgoingEdges(v).end(),
                             [&](EdgeId e) { return this->Strong(e); });
    }

    size_t DomainGraph::WeakEdgeCount(VertexId v) const {
        return std::count_if(this->OutgoingEdges(v).begin(), this->OutgoingEdges(v).end(),
                             [&](EdgeId e) { return !(this->Strong(e)); });
    }

    bool DomainGraph::HasStrongIncomingEdge(VertexId v) const {
        auto edges = this->IncomingEdges(v);
        for (auto e : edges) {
            if (this->Strong(e)) {
                return true;
            }
        }
        return false;
    }

    bool DomainGraph::NearContigStart(VertexId v) const {
        return this->data(v).GetNearContigStart();
    }

    bool DomainGraph::NearContigEnd(VertexId v) const {
        return this->data(v).GetNearContigEnd();
    }

    bool DomainGraph::Visited(VertexId v) const {
        return this->data(v).Visited();
    }

    std::vector<DomainGraph::EdgeId> DomainGraph::GetDomainEdges(VertexId v) const {
        return this->data(v).GetDomainEdges();
    }

    bool DomainGraph::Strong(EdgeId e) const {
        return this->data(e).Strong();
    }

    void DomainGraph::SetContigNearEnd(VertexId v) {
        this->data(v).SetNearEndCoord();
        this->data(this->conjugate(v)).SetNearStartCoord();
    }

    void DomainGraph::ExportToDot(const std::string &output_path) const {
        std::ofstream out(output_path);
        out << "digraph domain_graph {" << std::endl;
        for (auto iter = this->SmartVertexBegin(); ! iter.IsEnd(); ++iter) {
            VertexId v = (*iter);
            out << "\"" << this->GetVertexName(v) << "\""
                << " [label=\"" << this->data(v).GetType() << " " << this->GetVertexName(v) << " "
                << this->data(v).GetMaxVisited() << "\"];" << std::endl;
        }

        for (auto iter = this->SmartEdgeBegin(); ! iter.IsEnd(); ++iter) {
            EdgeId e = (*iter);
            out << "\"" << this->GetVertexName(this->EdgeStart(e)) << "\""
                << " -> "
                << "\"" << this->GetVertexName(this->EdgeEnd(e)) << "\"";
            out << " [label=";
            out << this->length(e);
            if (this->Strong(e)) {
                out << " style=bold];" << std::endl;
            } else {
                out << " style=dotted];" << std::endl;
            }
        }
        out << "}";
    }

    void DomainGraph::FindDomainOrderings(debruijn_graph::GraphPack &gp,
                             const std::string &output_filename, const std::string &output_dir) {
        const auto &graph = gp.get<debruijn_graph::Graph>();
        std::ofstream stat_stream(output_dir + "/bgc_statistics.txt");
        FindBasicStatistic(stat_stream);
        path_extend::ContigWriter writer(graph, path_extend::MakeContigNameGenerator(cfg::get().mode, gp));
        std::set<VertexId> nodes_with_incoming;
        for (auto iter = this->SmartEdgeBegin(); ! iter.IsEnd(); ++iter)
            nodes_with_incoming.insert(this->EdgeEnd(*iter));
        std::set<VertexId> start_nodes;
        for (auto iter = this->SmartVertexBegin(); ! iter.IsEnd(); ++iter) {
            VertexId v = *iter;
            if (nodes_with_incoming.find(v) == nodes_with_incoming.end() &&
                this->OutgoingEdgeCount(v)) {
                start_nodes.insert(v);
            }
        }
        io::osequencestream_bgc oss(output_dir + "/" + output_filename);
        path_extend::ScaffoldSequenceMaker seq_maker(graph);
        std::vector<std::vector<VertexId>> answer;
        int ordering_id = 1;
        int component_id = 1;

        for (auto v : start_nodes) {
            if (!this->Visited(v)) {
                stat_stream << "BGC subgraph " << component_id << std::endl;
                FindAllPossibleArrangements(v, answer, stat_stream);
                ordering_id = 1;
                for (auto vec : answer) {
                    OutputStatArrangement(vec, ordering_id, stat_stream);
                    path_extend::BidirectionalPath *p =
                            new path_extend::BidirectionalPath(graph);
                    std::string outputstring = PathToSequence(p, vec);
                    path_extend::BidirectionalPath *conjugate =
                            new path_extend::BidirectionalPath(p->Conjugate());
                    contig_paths_.AddPair(p, conjugate);
                    OutputComponent(p, component_id, ordering_id);
                    outputstring = seq_maker.MakeSequence(*p);
                    oss.SetCluster(component_id, ordering_id);
                    oss << outputstring;
                    ordering_id++;
                }
                answer.clear();
                component_id++;
            }
        }
    }
}
