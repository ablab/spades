//***************************************************************************
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/simplification/superbubble_finder.hpp"
#include "pipeline/graph_pack.hpp"
#include "visualization/visualization.hpp"
#include "visualization/graph_labeler.hpp"
#include "io/graph/gfa_writer.hpp"
#include "io/reads/io_helper.hpp"
#include "assembly_graph/graph_support/genomic_quality.hpp"

#include <unordered_set>

namespace bin_refinement {

using namespace debruijn_graph;

template<class SequenceMapper, class SingleReadStream>
double CountMedianCoverage(const Graph &g, const SequenceMapper &mapper, SingleReadStream &stream) {
    std::vector<std::pair<double, size_t>> cov_lens;
    size_t total = 0;
    io::SingleRead r;
    while (!stream.eof()) {
        stream >> r;
        for (auto e_mr : mapper.MapRead(r)) {
            cov_lens.push_back({g.coverage(e_mr.first), e_mr.second.mapped_range.size()});
            total += e_mr.second.mapped_range.size();
        }
    }
    INFO("Total length of contig mappings " << total);
    std::sort(cov_lens.begin(), cov_lens.end());
    size_t sum = 0;
    for (auto cov_len : cov_lens) {
        sum += cov_len.second;
        if (sum >= total / 2)
            return cov_len.first;
    }
    return 0.;
}

typedef std::set<EdgeId> EdgeSet;
typedef std::set<VertexId> VertexSet;

template<class NeighbourIteratorFactory>
class DFS {
    typedef vertex_neighbour<Graph> VNeighbour;
    typedef typename NeighbourIteratorFactory::NeighbourIterator NeighbourIterator;

    const Graph &g_;
    //traversal will start from the end of this edge
    const EdgeId init_edge_;
    const size_t edge_length_bound_;
    const size_t vertex_limit_;
    const VertexSet *blocked_ptr_;

    //todo switch to stack?
    std::vector<NeighbourIterator> edge_it_stack_;
    std::vector<VertexId> vertex_stack_;

    VertexSet reached_;
    EdgeSet border_sources_;
    EdgeSet border_sinks_;

    VertexId BeforeLastVertex() const {
        return vertex_stack_[vertex_stack_.size() - 2];
    }

public:

    DFS(const Graph &g, EdgeId e,
        size_t edge_length_bound,// = 10000,
        size_t vertex_limit = std::numeric_limits<size_t>::max(),
        const VertexSet *blocked_ptr = nullptr) :
            g_(g), init_edge_(e),
            edge_length_bound_(edge_length_bound),
            vertex_limit_(vertex_limit),
            blocked_ptr_(blocked_ptr) {

        vertex_stack_.push_back(g_.EdgeStart(e));
        VertexId v = g_.EdgeEnd(e);
        DEBUG("Initializing: edge " << g_.int_id(e) << " vertex " << g_.int_id(v));
        if (!blocked_ptr_ || !blocked_ptr_->count(v)) {
            reached_.insert(v);
            DEBUG("Initializing: pushing neighbour iterator for vertex " << g_.int_id(v));
            edge_it_stack_.emplace_back(NeighbourIterator(g_, v));
            vertex_stack_.push_back(v);
            DEBUG("Border source edge");
            border_sources_.insert(e);
        }
    }

    //returns false if vertex limit was exceeded
    //FIXME should long loop (etc) within component be source/sink?
    bool Run() {
        while (!edge_it_stack_.empty()) {
            auto &it = edge_it_stack_.back();
            if (it.HasNext()) {
                VNeighbour neighbour = it.Next();
                VertexId v = neighbour.vertex;
                EdgeId e = neighbour.edge;
                DEBUG("Considering neighbour edge " << g_.str(e) << " vertex " << g_.int_id(v));

                if (g_.IncomingEdgeCount(v) == 0) {
                    DEBUG("Border source edge");
                    border_sources_.insert(e);
                    continue;
                }
                if (g_.OutgoingEdgeCount(v) == 0) {
                    DEBUG("Border sink edge");
                    border_sinks_.insert(e);
                    continue;
                }

                if (g_.length(e) >= edge_length_bound_) {
                    if (g_.EdgeEnd(e) == v) {
                        DEBUG("Border sink edge");
                        border_sinks_.insert(e);
                    }
                    if (g_.EdgeStart(e) == v) {
                        DEBUG("Border source edge");
                        border_sources_.insert(e);
                    }

                    continue;
                }

                ////just for premature optimization reasons :)
                //if (v == BeforeLastVertex()) {
                //    //not going right back from we came from
                //    DEBUG("Preventive optimization for bidirectional case");
                //    continue;
                //}

                if (blocked_ptr_ && blocked_ptr_->count(v)) {
                    DEBUG("Blocked");
                    continue;
                }

                if (!reached_.count(v)) {
                    //fixme Process v
                    reached_.insert(v);
                    if (reached_.size() > vertex_limit_) {
                        DEBUG("Vertex limit reached");
                        return false;
                    }
                    DEBUG("Pushing neighbour iterator for vertex " << g_.int_id(v));
                    edge_it_stack_.emplace_back(NeighbourIterator(g_, v));
                    //fixme remove vertex_stack_?
                    vertex_stack_.push_back(v);
                } else {
                    DEBUG("Already considered");
                }
            } else {
                edge_it_stack_.pop_back();
                vertex_stack_.pop_back();
            }
        }
        return true;
    }

    const EdgeSet &border_sources() const {
        return border_sources_;
    }

    const EdgeSet &border_sinks() const {
        return border_sinks_;
    }

    const VertexSet &reached() const {
        return reached_;
    }

};

/* Expands on edges located between long annotated ones.
 * Let's call path consisting of "short" edges an s.e.p.
 * In particular we add all edges induced by vertices, s.t.:
 * 1. v can be accessed via an s.e.p. from long annotated edge
 * 2. v can not be accessed via an s.e.p. from long non-annotated edge
 * 3. there is an s.e.p. connecting v to a long annotated edge
 * 4. there is no s.e.p. connecting v to a long non-annotated edge
 */
//FIXME add size check
//TODO what about edges < 10Kbp annotated with a different bin?
class BoundedRegionsExpander {
    const Graph &g_;
    const std::set<EdgeId> &annotated_;
    const size_t edge_length_bound_;

    bool LongEdge(EdgeId e) const {
        return g_.length(e) >= edge_length_bound_;
    }

    EdgeSet RCEdges(const EdgeSet &edges) const {
        EdgeSet answer;
        for (EdgeId e : edges) {
            answer.insert(g_.conjugate(e));
        }
        return answer;
    }

    VertexSet RCEdges(const VertexSet &vertices) const {
        VertexSet answer;
        for (VertexId v : vertices) {
            answer.insert(g_.conjugate(v));
        }
        return answer;
    }

    VertexSet CollectBadVerticesFWD(const EdgeSet &border_sources) const {
        typedef DFS<ForwardNeighbourIteratorFactory<Graph>> ForwardDFS;
        VertexSet bad;
        for (EdgeId e : border_sources) {
            if (!annotated_.count(e)) {
                VERIFY(!annotated_.count(g_.conjugate(e)));
                ForwardDFS dfs(g_, e, edge_length_bound_,
                        std::numeric_limits<size_t>::max(), &bad);
                dfs.Run();
                utils::insert_all(bad, dfs.reached());
            }
        }
        return bad;
    }

    VertexSet CollectGoodVerticesFWD(const EdgeSet &border_sources, const VertexSet &bad) const {
        typedef DFS<ForwardNeighbourIteratorFactory<Graph>> ForwardDFS;
        VertexSet good;
        for (EdgeId e : border_sources) {
            if (annotated_.count(e)) {
                VERIFY(annotated_.count(g_.conjugate(e)));
                ForwardDFS dfs(g_, e, edge_length_bound_,
                        std::numeric_limits<size_t>::max(), &bad);
                dfs.Run();
                utils::insert_all(good, dfs.reached());
            }
        }
        return good;
    }

    //FIXME can be simplified with GraphComponent?
    //FIXME is it ok that the set of reached vertices was not used?!
    //It also should find sources and sinks automatically?!
    EdgeSet ProcessComponent(const VertexSet &reached,
                             const EdgeSet &border_sources,
                             const EdgeSet &border_sinks) const {
        //FIXME can I optimize further by collecting fwd/bwd bad vertices first?
        VertexSet good_from_sources = CollectGoodVerticesFWD(border_sources, CollectBadVerticesFWD(border_sources));
        auto rc_sinks = RCEdges(border_sinks);
        VertexSet good_from_rc_sinks = CollectGoodVerticesFWD(rc_sinks, CollectBadVerticesFWD(rc_sinks));
        EdgeSet to_expand;
        for (VertexId v : good_from_sources) {
            for (EdgeId e : g_.OutgoingEdges(v)) {
                if (good_from_rc_sinks.count(g_.conjugate(g_.EdgeEnd(e)))) {
                    VERIFY(reached.count(v) && reached.count(g_.EdgeEnd(e)));
                    to_expand.insert(e);
                }
            }
        }
        return to_expand;
    }

public:

    //FIXME check that annotated set contains rc edges
    BoundedRegionsExpander(const Graph &g,
                           const std::set<EdgeId> &annotated,
                           size_t edge_length_bound/* = 10000*/) :
            g_(g), annotated_(annotated),
            edge_length_bound_(edge_length_bound) {
        for (EdgeId e : annotated) {
            VERIFY(annotated.count(g.conjugate(e)) > 0);
        }
    }

    //FIXME don't forget to exclude vertices of annotated border edges from set of bad vertices
    std::set<EdgeId> Run() const {
        typedef DFS<UnorientedNeighbourIteratorFactory<Graph>> ComponentCollectingDFS;
        std::set<EdgeId> considered;
        std::set<EdgeId> to_expand;
        for (EdgeId e : annotated_)  {
            DEBUG("Analyzing region for the end of edge " << g_.str(e));
            if (LongEdge(e) && !considered.count(e)) {
                //FIXME length bounds have to by synched here!!!
                ComponentCollectingDFS component_collector(g_, e, edge_length_bound_);
                component_collector.Run();

                //FIXME "reached" not used in any meaningful way
                DEBUG("Considering short-edge component of size " << component_collector.reached().size());
                utils::insert_all(to_expand,
                                  ProcessComponent(component_collector.reached(),
                                                   component_collector.border_sources(),
                                                   component_collector.border_sinks()));

                utils::insert_all(considered, component_collector.border_sources());
                utils::insert_all(considered, RCEdges(component_collector.border_sinks()));
            }
        }
        return to_expand;
    }
};

//static size_t TotalLength(const Graph &g, const EdgeSet & edges) {
//    size_t tot_len = 0;
//    for (EdgeId e : edges) {
//        tot_len += g.length(e);
//    }
//    return tot_len;
//}

static EdgeSet ExtraEdges(const Graph &g, const EdgeSet &expanded, const EdgeSet &annotated) {
    EdgeSet extra;
    for (EdgeId e : expanded) {
        if (!annotated.count(e) && !extra.count(e)) {
            VERIFY(!annotated.count(g.conjugate(e)));
            //DEBUG("Adding edge to the annotation expansion set" << g.str(e));
            extra.insert(e);
            extra.insert(g.conjugate(e));
        }
    }
    return extra;
}

static void DrawComponent(const Graph &g, const VertexSet &component_vertices, const std::string &path,
        const EdgeSet &annotated, const EdgeSet &to_expand = EdgeSet(),
        const omnigraph::EdgesPositionHandler<Graph> *position_handler = nullptr) {
    VertexSet neighbourhood(component_vertices);
    for (VertexId v : component_vertices) {
        for (EdgeId edge : g.IncidentEdges(v)) {
            neighbourhood.insert(g.EdgeStart(edge));
            neighbourhood.insert(g.EdgeEnd(edge));
        }
    }
    auto gc = GraphComponent<Graph>::FromVertices(g, neighbourhood);
//    DEBUG("Writing to gfa");
//    std::ofstream os(pics_path_ + std::to_string(g_.int_id(e)) + ".gfa");
//    gfa::GFAWriter writer(g_, os);
//    writer.WriteSegmentsAndLinks(gc);
    using namespace visualization::graph_colorer;
    auto edge_colorer = std::make_shared<CompositeEdgeColorer<Graph>>("black");
    edge_colorer->AddColorer(std::make_shared<SetColorer<Graph>>(g, annotated, "blue"));
    edge_colorer->AddColorer(std::make_shared<SetColorer<Graph>>(g, to_expand, "red"));
    visualization::visualization_utils::WriteComponent(gc,
                                                       path,
            //visualization::graph_colorer::DefaultColorer(g_),
                                                       DefaultColorer(g, edge_colorer),
                                                       visualization::graph_labeler::DefaultLabeler<Graph>(g, position_handler));
}

//FIXME do not add long edges of the induced subgraph?!!

/* Expands on edges for which all the long edges it is reachable from are annotated.
 * Let's call path consisting of "short" edges an s.e.p.
 * In particular we add all edges, s.t. for their starting vertex v:
 * 1. v can be accessed via an s.e.p. from long annotated edge
 * 2. v can not be accessed via an s.e.p. from long non-annotated edge
 */
class OnlyAnnotatedReachableExpander {
    const Graph &g_;
    const EdgeSet &annotated_;
    const size_t edge_length_bound_;

    //TODO think if we need it. Currently always equal to edge_length_bound
    const size_t coverage_check_length_bound_;
    const double cov_lower_bound_;
    const double cov_upper_bound_;
    const std::string &pics_path_;

    bool LongEdge(EdgeId e) const {
        return g_.length(e) >= edge_length_bound_;
    }

    EdgeSet RCEdges(const EdgeSet &edges) const {
        EdgeSet answer;
        for (EdgeId e : edges) {
            answer.insert(g_.conjugate(e));
        }
        return answer;
    }

    VertexSet RCEdges(const VertexSet &vertices) const {
        VertexSet answer;
        for (VertexId v : vertices) {
            answer.insert(g_.conjugate(v));
        }
        return answer;
    }

    //Collecting vertices reachable from unannotated border sources as "bad" vertices
    VertexSet CollectBadVerticesFWD(const EdgeSet &border_sources) const {
        typedef DFS<ForwardNeighbourIteratorFactory<Graph>> ForwardDFS;
        VertexSet bad;
        for (EdgeId e : border_sources) {
            if (!annotated_.count(e)) {
                INFO("Unannotated border edge " << g_.str(e));
                VERIFY(!annotated_.count(g_.conjugate(e)));
                ForwardDFS dfs(g_, e, edge_length_bound_, std::numeric_limits<size_t>::max(), &bad);
                dfs.Run();
                utils::insert_all(bad, dfs.reached());
            }
        }
        return bad;
    }

    //Collecting vertices reachable from annotated border sources.
    //Blocking on "bad" vertices
    VertexSet CollectGoodVerticesFWD(const EdgeSet &border_sources, const VertexSet &bad) const {
        typedef DFS<ForwardNeighbourIteratorFactory<Graph>> ForwardDFS;
        VertexSet good;
        for (EdgeId e : border_sources) {
            if (annotated_.count(e)) {
                INFO("Annotated border edge " << g_.str(e));
                VERIFY(annotated_.count(g_.conjugate(e)));
                ForwardDFS dfs(g_, e, edge_length_bound_, std::numeric_limits<size_t>::max(), &bad);
                dfs.Run();
                utils::insert_all(good, dfs.reached());
            }
        }
        return good;
    }

    //FIXME is it ok that the set of reached vertices was not used?!
    EdgeSet ProcessComponent(const VertexSet &/*reached*/,
                             const EdgeSet &border_sources) const {
        VertexSet good_from_sources = CollectGoodVerticesFWD(border_sources, CollectBadVerticesFWD(border_sources));
        INFO("Considering edges from " << good_from_sources.size() << " vertices");
        EdgeSet to_expand;
        for (VertexId v : good_from_sources) {
            //TODO improve checks.
            // E.g. do not go beyond edges for which check failed.
            for (EdgeId e : g_.OutgoingEdges(v)) {
                if (g_.length(e) >= coverage_check_length_bound_) {
                    if (math::ls(g_.coverage(e), cov_lower_bound_)
                        || math::gr(g_.coverage(e), cov_upper_bound_)) {
                        continue;
                    }
                }
                to_expand.insert(e);
            }
        }
        return to_expand;
    }

public:

    //FIXME check that annotated set contains rc edges
    //NB: cov_bounds should be disabled while working with strain mixture graphs
    OnlyAnnotatedReachableExpander(const Graph &g,
                                   const EdgeSet &annotated,
                                   size_t edge_length_bound/* = 10000*/,
                                   double cov_lower_bound,
                                   double cov_upper_bound,
                                   const std::string &pics_folder = "") :
            g_(g), annotated_(annotated),
            edge_length_bound_(edge_length_bound),
            coverage_check_length_bound_(edge_length_bound),
            cov_lower_bound_(cov_lower_bound),
            cov_upper_bound_(cov_upper_bound),
            pics_path_(pics_folder) {
        INFO("Edge length bound " << edge_length_bound);
        for (EdgeId e : annotated) {
            VERIFY(annotated.count(g.conjugate(e)) > 0);
        }
    }

    //TODO don't forget to exclude vertices of annotated border edges from set of bad vertices
    //Don't understand this comment any more :)
    std::set<EdgeId> Run() const {
        typedef DFS<UnorientedNeighbourIteratorFactory<Graph>> ComponentCollectingDFS;
        std::set<EdgeId> considered;
        std::set<EdgeId> answer;
        for (EdgeId e : annotated_) {
            if (LongEdge(e) && !considered.count(e)) {
                DEBUG("Processing edge " << g_.str(e));

                ComponentCollectingDFS component_collector(g_, e, edge_length_bound_, /*vertex limit*/500);
                bool limit_not_exceeded = component_collector.Run();

                utils::insert_all(considered, component_collector.border_sources());

                if (!limit_not_exceeded) {
                    DEBUG("Component collected from edge " << g_.str(e) << " was too large");
                    continue;
                }

                //Using less straightforward code for visualization
//                utils::insert_all(answer,
//                                  ProcessComponent(component_collector.reached(),
//                                                   component_collector.border_sources()));

                auto extra = ExtraEdges(g_,
                                            ProcessComponent(component_collector.reached(),
                                                    component_collector.border_sources()),
                                            annotated_);
                DEBUG("Size of set to expand " << extra.size());

                if (!pics_path_.empty() && extra.size() > 0) {
                    DEBUG("Writing to dot");
                    DrawComponent(g_,
                            component_collector.reached(), pics_path_ + std::to_string(g_.int_id(e)) + ".dot",
                            annotated_, extra);
                    DEBUG("Written");
                }


                DEBUG("Done processing edge " << g_.str(e));
                utils::insert_all(answer, extra);
            }
        }
        return answer;
    }

private:
    DECL_LOGGER("OnlyAnnotatedReachableExpander");
};

static EdgeSet UnambiguousExpand(const Graph &g, VertexId v) {
    EdgeSet expanded;
    DEBUG("Unambiguously extending vertex " << g.str(v));
    while (true) {
        DEBUG("Checking vertex " << g.str(v));
        if (g.OutgoingEdgeCount(v) == 1) {
            DEBUG("Single outgoing edge");
            EdgeId e = g.GetUniqueOutgoingEdge(v);
            if (expanded.count(e) > 0) {
                //return if already annotated
                //prevents looping and processing same edges multiple times
                return expanded;
            }
            expanded.insert(e);
            v = g.EdgeEnd(e);
            continue;
        }

        omnigraph::SuperbubbleFinder<Graph> superbubble_finder(g, v);
        DEBUG("Superbubble search");
        if (superbubble_finder.FindSuperbubble()) {
            auto gc = superbubble_finder.AsGraphComponent();
            DEBUG("Superbubble: v_size " << gc.v_size() << " e_size " << gc.e_size());

            if (std::all_of(gc.e_begin(), gc.e_end(), [&] (EdgeId e) {
                        return expanded.count(e) > 0;})) {
                //exit if the entire component already annotated
                //prevents looping and processing same edges multiple times
                DEBUG("Obsolete bubble");
                return expanded;
            }

            DEBUG("Adding " << gc.edges().size() << " edges");
            for (EdgeId e : gc.edges()) {
                expanded.insert(e);
            }
            v = superbubble_finder.end_vertex();
            continue;
        }
        return expanded;
    }
}

//Fixme introduce length condition?
static EdgeSet UnambiguousExtensions(const Graph &g,
                                              const EdgeSet &annotated,
                                              const std::string &pics_path = "") {
    INFO("Unambiguous extension started");
    std::set<EdgeId> all_extra;
    for (EdgeId e : annotated) {
        VertexId v = g.EdgeEnd(e);
        EdgeSet expanded = UnambiguousExpand(g, v);
        auto extra = ExtraEdges(g, expanded, annotated);

        if (!pics_path.empty() && extra.size() > 0) {
            DEBUG("Writing to dot");
            DrawComponent(g,
                          GraphComponent<Graph>::FromEdges(g, expanded.begin(), expanded.end()).vertices(),
                          pics_path + std::to_string(g.int_id(v)) + ".dot",
                          annotated, extra);
            DEBUG("Written");
        }
        utils::insert_all(all_extra, extra);
    }
    INFO("Unambiguous extension done");
    return all_extra;
}

static void AnnotateExtraEdgesRound(const Graph &g,
                                    EdgeQuality<Graph> &edge_qual,
                                    double est_cov,
                                    const std::string &out_path) {
    INFO("Annotation expansion");
    auto annotated_edges = edge_qual.PositiveQualEdges();

    fs::make_dirs(out_path + "/unambig/");
    fs::make_dirs(out_path + "/only_ann_reach/");
    EdgeSet expanded(annotated_edges);
    utils::insert_all(expanded, UnambiguousExtensions(g, annotated_edges, out_path + "/unambig/"));
    ////TODO think about step by step additions (currently initial annotation set is not changed)
    //utils::insert_all(expanded, BoundedRegionsExpander(g, annotated_edges, 10000).Run());
    //INFO("Extra edges total " << expanded.size());
    utils::insert_all(expanded, OnlyAnnotatedReachableExpander(g, annotated_edges,
            //FIXME magic constants
                                                                  2500, 0.7 * est_cov, 1.3 * est_cov,
                                                                  out_path + "/only_ann_reach/").Run());
    INFO("Only annotated reachable done");

    std::set<EdgeId> extra_edges;
    size_t extra_cnt = 0;
    size_t total_length = 0;
    for (EdgeId e : expanded) {
        if (!annotated_edges.count(e) && !extra_edges.count(e)) {
            VERIFY(!annotated_edges.count(g.conjugate(e)));
            extra_cnt++;
            total_length += g.length(e);
            extra_edges.insert(e);
            extra_edges.insert(g.conjugate(e));
        }
    }
    INFO("Extra edges: " << extra_cnt << "; of length " << total_length);

    static const double QUALITY_MARKER = 1000.;
    for (EdgeId e : extra_edges) {
        edge_qual.AddQuality(e, QUALITY_MARKER);
    }
    INFO("Annotation expansion done");
}

static void CompareAnnotations(const GraphPack &gp, const EdgeQuality<Graph> &annotation,
                        const std::string &out_folder = "") {
    const auto &g = gp.get<Graph>();
    const auto &edge_qual = gp.get<EdgeQuality<Graph>>();
    using namespace visualization::graph_colorer;
    using namespace visualization::graph_labeler;
    using namespace visualization::visualization_utils;
    auto edge_colorer = std::make_shared<CompositeEdgeColorer<Graph>>("black");
    if (!out_folder.empty()) {
        edge_colorer->AddColorer(std::make_shared<SetColorer<Graph>>(g, edge_qual.PositiveQualEdges(), "blue"));
        edge_colorer->AddColorer(std::make_shared<SetColorer<Graph>>(g, annotation.PositiveQualEdges(), "red"));
    }

    const size_t EDGE_LEN = 10000;

    const auto &edge_pos= gp.get<EdgesPositionHandler<Graph>>();
    size_t total_len = 0;
    for (EdgeId e : edge_qual.PositiveQualEdges()) {
        if (g.int_id(e) <= g.int_id(g.conjugate(e)) && annotation.IsZeroQuality(e)) {
            if (g.length(e) > EDGE_LEN) {
                INFO("Reference edge " << g.str(e) << " missing from the annotation");
                if (!out_folder.empty()) {
                    DefaultLabeler<Graph> labeler(g, edge_pos);
                    edge_colorer->AddColorer(std::make_shared<SetColorer<Graph>>(g, std::vector<EdgeId>(1, e), "green"));

                    WriteComponent<Graph>(omnigraph::VertexNeighborhood<Graph>(g, g.EdgeStart(e), 50, 5000),
                                          out_folder + "/edge_" + std::to_string(g.int_id(e)) + "_s.dot",
                                          DefaultColorer(g, edge_colorer), labeler);
                    WriteComponent<Graph>(omnigraph::VertexNeighborhood<Graph>(g, g.EdgeEnd(e), 50, 5000),
                                          out_folder + "/edge_" + std::to_string(g.int_id(e)) + "_e.dot",
                                          DefaultColorer(g, edge_colorer), labeler);
                    edge_colorer->PopColorer();
                }
            }
            total_len += g.length(e);
        }
    }
    INFO("Total length " << total_len);
}

static void AnnotateExtraEdges(const GraphPack &gp,
                               EdgeQuality<Graph> &annotation,
                               double est_cov,
                               const std::string &out_folder = "") {
    const size_t expansion_it_cnt = 3;
    const auto &edge_qual = gp.get<EdgeQuality<Graph>>();
    const auto &graph = gp.get<Graph>();

    for (size_t i = 0 ; i < expansion_it_cnt ; ++i) {
        if (edge_qual.IsAttached()) {
            INFO("Reference was provided. Comparing annotations.");

            CompareAnnotations(gp, annotation);
        }

        INFO("#" << (i + 1) << " round of expansion")
        AnnotateExtraEdgesRound(graph, annotation, est_cov, out_folder);

    }

    if (edge_qual.IsAttached()) {
        INFO("Reference was provided. Comparing annotations.");
        fs::make_dirs(out_folder + "/unannotated_pics");
        CompareAnnotations(gp, annotation, out_folder + "/unannotated_pics");
    }
}

static double AnnotateEdges(const GraphPack &gp,
                           EdgeQuality<Graph> &annotation,
                           const std::string &bin_contigs) {
    auto stream = io::EasyStream(bin_contigs,
            /*followed by rc*/true);
    const auto &graph = gp.get<Graph>();
    const auto &index = gp.get<EdgeIndex<Graph>>();
    const auto &kmer_mapper = gp.get<KmerMapper<Graph>>();
    BasicSequenceMapper<Graph, EdgeIndex<Graph>> mapper(graph, index, kmer_mapper,
                                                        /*extension optimization enabled*/false);

    INFO("Filling annotation with contigs attributed to bin");
    annotation.Fill(mapper, stream);

    VERIFY(annotation.IsAttached());
    stream.reset();

    double base_cov = CountMedianCoverage(graph, mapper, stream);

    INFO("Base coverage determined as " << base_cov);
    return base_cov;
}

}
