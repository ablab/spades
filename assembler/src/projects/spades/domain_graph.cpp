//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "domain_graph.hpp"

#include "io/reads/osequencestream.hpp"
#include "pipeline/graph_pack.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <memory>

using namespace nrps;

void DomainGraph::ExportToDot(const std::string &output_path) const {
    std::ofstream out(output_path);
    out << "digraph domain_graph {" << std::endl;
    for (auto v : nodes_) {
        out << "\"" << v->name_ << "\""
            << " [label=\"" << v->v_type_ << " " << v->name_ << " "
            << v->max_visited_ << "\"];" << std::endl;
    }

    for (auto e : arcs_) {
        out << "\"" << e->start_->name_ << "\""
            << " -> "
            << "\"" << e->end_->name_ << "\"";
        out << " [label=";
        out << e->length_;
        if (e->strong_) {
            out << " style=bold];" << std::endl;
        } else {
            out << " style=dotted];" << std::endl;
        }
    }
    out << "}";
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

    for (auto v : nodes_)
        domain_count[v->v_type_]++;

    for (auto domain_type : domain_count) {
        stat_stream << "# " << domain_type.first << "-domains - "
                    << domain_type.second << std::endl;
    }
    stat_stream << std::endl;
}

void DomainGraph::FindDomainOrderings(debruijn_graph::conj_graph_pack &gp,
                                      const std::string &output_filename, const std::string &output_dir) {
    std::ofstream stat_stream(output_dir + "/bgc_statistics.txt");
    FindBasicStatistic(stat_stream);
    path_extend::ContigWriter writer(gp.g, path_extend::MakeContigNameGenerator(cfg::get().mode, gp));
    std::set<std::shared_ptr<Vertex>> nodes_with_incoming;
    for (auto e : arcs_)
        nodes_with_incoming.insert(e->end_);

    std::set<std::shared_ptr<Vertex>> start_nodes;
    for (auto v : nodes_) {
        if (nodes_with_incoming.find(v) == nodes_with_incoming.end() &&
            v->edges_.size()) {
            start_nodes.insert(v);
        }
    }
    io::osequencestream_bgc oss(output_dir + "/" + output_filename);
    path_extend::ScaffoldSequenceMaker seq_maker(gp.g);
    std::vector<std::vector<std::shared_ptr<Vertex>>> answer;
    int ordering_id = 1;
    int component_id = 1;

    for (auto v : start_nodes) {
        if (!v->visited_) {
            stat_stream << "BGC subgraph " << component_id << std::endl;
            FindAllPossibleArrangements(gp.g, v, answer, stat_stream);
            ordering_id = 1;
            for (auto vec : answer) {
                OutputStatArrangement(gp.g, vec, ordering_id, stat_stream);
                path_extend::BidirectionalPath *p =
                        new path_extend::BidirectionalPath(gp.g);
                std::string outputstring = PathToSequence(p, vec);
                path_extend::BidirectionalPath *conjugate =
                        new path_extend::BidirectionalPath(p->Conjugate());
                contig_paths_.AddPair(p, conjugate);
                OutputComponent(gp, p, component_id, ordering_id);
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

void DomainGraph::ExportPaths(debruijn_graph::conj_graph_pack &gp,
                              const std::string &output_dir) {
    auto name_generator =
            path_extend::MakeContigNameGenerator(cfg::get().mode, gp);
    path_extend::ContigWriter writer(gp.g, name_generator);
    writer.OutputPaths(contig_paths_, output_dir + "/orderings2");
}

static std::set<EdgeId> CollectEdges(path_extend::BidirectionalPath *p) {
    std::set<EdgeId> edge_set;
    for (size_t i = 0; i < p->Size(); ++i) {
        edge_set.insert(p->At(i));
    }
    return edge_set;
}

void DomainGraph::OutputComponent(debruijn_graph::conj_graph_pack &gp,
                                  path_extend::BidirectionalPath *p, int component_id,
                                  int ordering_id) {
    auto edges = CollectEdges(p);
    auto comp = GraphComponent<Graph>::FromEdges(gp.g, edges.begin(),
                                                 edges.end(), true);
    std::ofstream os(cfg::get().output_dir + "/bgc_in_gfa/" +
                     std::to_string(component_id) + "_" + std::to_string(ordering_id) +
                     ".gfa");
    path_extend::GFAWriter<Graph> writer(gp.g, os);
    writer.Write(comp);
}

std::string DomainGraph::PathToSequence(path_extend::BidirectionalPath *p,
                                        std::vector<std::shared_ptr<Vertex>> &answer) {
    std::stringstream ss;
    DEBUG("Translating " << p->GetId() << " to sequnece");
    for (size_t i = 0; i < answer.size(); ++i) {

        auto v = answer[i];
        DEBUG("Translating vertex " << v->name_);
        
        for (auto e : v->domain_edges_in_row_) {
            if (p->Size() == 0 || p->Back() != e) {
                int gap = 0;
                if (p->Size() != 0 &&
                    p->graph().EdgeEnd(p->Back()) != p->graph().EdgeStart(e)) {
                    gap = 100;
                }
                p->PushBack(e, path_extend::Gap(gap));
            }
        }
        if (i != answer.size() - 1) {
            Edge *next_edge = this->getEdge(v, answer[i + 1]);
            DEBUG("Translating edge between " <<
                  v->name_ << " and " << answer[i + 1]->name_);
            DEBUG("Is edge strong " << next_edge->strong_);
            if (next_edge == nullptr) {
                DEBUG("Something strange!");
                continue;
            }
            for (auto e : next_edge->edges_) {
                if (p->Size() == 0 || p->Back() != e) {
                    int gap = 0;
                    if (p->Size() != 0 &&
                        p->graph().EdgeEnd(p->Back()) != p->graph().EdgeStart(e)) {
                        gap = 100;
                    }
                    p->PushBack(e, path_extend::Gap(gap));
                }
            }
        }
    }
    return ss.str();
}

void DomainGraph::OutputStatArrangement(const Graph &g,
                                        std::vector<std::shared_ptr<Vertex>> single_candidate,
                                        int id, std::ofstream &stat_file) {
    stat_file << "BGC candidate " << id << std::endl;
    std::string delimeter = "";
    bool is_nrps = false;
    bool is_pks = false;
    for (auto v : single_candidate) {
        stat_file << delimeter;
        delimeter = "-";
        if (v->v_type_ == "AMP") {
            stat_file << "A";
            is_nrps = true;
        } else if (v->v_type_ == "CStart") {
            stat_file << "C";
            is_nrps = true;
        } else if (v->v_type_ == "AT") {
            stat_file << "AT";
            is_pks = true;
        } else if (v->v_type_ == "TE") {
            stat_file << "TE";
        } else if (v->v_type_ == "KR") {
            stat_file << "KR";
            is_pks = true;
        } else if (v->v_type_ == "KS") {
            stat_file << "KS";
            is_pks = true;
        } else {
            stat_file << v->v_type_;
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
        stat_file << "None";
    }
    stat_file << std::endl;
    stat_file << "Domain cordinates:" << std::endl;
    size_t current_coord = 0;
    path_extend::BidirectionalPath *p = new path_extend::BidirectionalPath(g);
    for (size_t i = 0; i < single_candidate.size(); ++i) {
        auto v = single_candidate[i];
        DEBUG("Translating vertex " << v->name_);
        for (auto e : v->domain_edges_in_row_) {
            if (p->Size() == 0 || p->Back() != e) {
                int gap = 0;
                if (p->Size() != 0 &&
                    p->graph().EdgeEnd(p->Back()) != p->graph().EdgeStart(e)) {
                    gap = 100;
                }
                p->PushBack(e, path_extend::Gap(gap));
                current_coord += g.length(e) + gap;
            }
        }

        size_t sum = 0;
        for (auto e2 : v->domain_edges_in_row_) {
            sum += g.length(e2);
        }

        stat_file << v->start_coord_ + current_coord - sum << " ";

        stat_file << v->end_coord_ + current_coord - g.length(p->Back())
                  << std::endl;

        if (i != single_candidate.size() - 1) {
            Edge *next_edge = this->getEdge(v, single_candidate[i + 1]);
            if (next_edge == nullptr) {
                continue;
            }
            for (auto e : next_edge->edges_) {
                if (p->Size() == 0 || p->Back() != e) {
                    int gap = 0;
                    if (p->Size() != 0 &&
                        p->graph().EdgeEnd(p->Back()) != p->graph().EdgeStart(e)) {
                        gap = 100;
                    }
                    p->PushBack(e, path_extend::Gap(gap));
                    current_coord += g.length(e) + gap;
                }
            }
        }
    }
}

void DomainGraph::OutputStat(std::set<std::shared_ptr<Vertex>> &preliminary_visited,
                             std::ofstream &stat_file) {
    stat_file << "# domains in the component - ";
    stat_file << preliminary_visited.size() << std::endl;
    int strong_edge_count =
            std::accumulate(preliminary_visited.begin(), preliminary_visited.end(),
                            0, [&](int prev, std::shared_ptr<Vertex> v) {
                                return prev + StrongEdgeCount(v);
                            });
    int weak_edge_count =
            std::accumulate(preliminary_visited.begin(), preliminary_visited.end(),
                            0, [&](int prev, std::shared_ptr<Vertex> v) {
                                return prev + (WeakEdgeCount(v));
                            });

    stat_file << "# strong/weak edges in the component - ";
    stat_file << strong_edge_count << "/" << weak_edge_count << std::endl;

    if (preliminary_visited.size() <= 5) {
        //                stat_file << "Warning: small number of domains, possible
        //                fragmentation" << std::endl;
    }

    if (strong_edge_count < (int)preliminary_visited.size() / 2) {
        //                stat_file << "Warning: small number of strong edges,
        //                unlikely to recover complete BGC" << std::endl;
    }
}

void DomainGraph::FindAllPossibleArrangements(const Graph &g, std::shared_ptr<Vertex> v,
                                              std::vector<std::vector<std::shared_ptr<Vertex>>> &answer,
                                              std::ofstream &stat_file) {
    std::set<std::shared_ptr<Vertex>> preliminary_visited;
    preliminary_visited.insert(v);
    PrelimDFS(v, preliminary_visited);
    size_t component_size = preliminary_visited.size();
    if (component_size == 1)
        return;

    SetCopynumber(g, preliminary_visited);
    OutputStat(preliminary_visited, stat_file);
    std::vector<std::shared_ptr<Vertex>> current;
    for (auto v : preliminary_visited) {
        v->visited_ = true;
        v->rc_->visited_ = true;
    }
    preliminary_visited.clear();
    iteration_number_ = 0;
    FinalDFS(v, current, preliminary_visited, answer, component_size);
}

void DomainGraph::PrelimDFS(std::shared_ptr<Vertex> v,
                            std::set<std::shared_ptr<Vertex>> &preliminary_visited) {
    for (auto e : v->edges_) {
        std::shared_ptr<Vertex> to = e->end_;
        if (preliminary_visited.find(to) == preliminary_visited.end()) {
            preliminary_visited.insert(to);
            PrelimDFS(to, preliminary_visited);
        }
    }
}

void DomainGraph::FinalDFS(std::shared_ptr<Vertex> v,
                           std::vector<std::shared_ptr<Vertex>> &current,
                           std::set<std::shared_ptr<Vertex>> preliminary_visited,
                           std::vector<std::vector<std::shared_ptr<Vertex>>> &answer,
                           size_t component_size) {
    iteration_number_++;

    current.push_back(v);
    v->visited_times_++;
    if (answer.size() > 50 || iteration_number_ > 10000000)
        return;

    bool was_extended = false;
    if (this->HasStrongEdge(v)) {
        for (auto arc : v->edges_) {
            if (arc->strong_) {
                std::shared_ptr<Vertex> to = arc->end_;
                if (preliminary_visited.find(to) != preliminary_visited.end())
                    continue;

                was_extended = true;
                FinalDFS(to, current, preliminary_visited, answer, component_size);
            }
        }
    } else {
        for (auto arc : v->edges_) {
            std::shared_ptr<Vertex> to = arc->end_;
            if (to->visited_times_ >= to->max_visited_)
                continue;

            if (preliminary_visited.find(to) != preliminary_visited.end())
                continue;

            was_extended = true;
            FinalDFS(to, current, preliminary_visited, answer, component_size);
        }
    }
    if (current.size() >= component_size && was_extended == false) {
        answer.push_back(current);
    }
    v->visited_times_--;
    current.pop_back();
}
