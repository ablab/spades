//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/simple_tools.hpp"

//#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
//#include "io/reads_io/file_reader.hpp"
#include "modules/simplification/tip_clipper.hpp"
#include "propagate.hpp"
#include "visualization.hpp"

namespace debruijn_graph {
static const size_t EDGE_LENGTH_THRESHOLD = 2000;

//FIXME 2kb edge length threshold might affect tip propagator in undesired way
class EdgeAnnotationPropagator {
    const conj_graph_pack& gp_;
    const string name_;
    size_t edge_length_threshold_;

protected:
    const conj_graph_pack& gp() const {
        return gp_;
    }

    const Graph& g() const {
        return gp_.g;
    }

    virtual set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const = 0;

public:
    EdgeAnnotationPropagator(const conj_graph_pack& gp,
                             const string& name,
                             size_t edge_length_threshold = EDGE_LENGTH_THRESHOLD) :
                    gp_(gp),
                    name_(name),
                    edge_length_threshold_(edge_length_threshold) {}

    const std::string& name() const {
        return name_;
    }

    std::map<bin_id, set<EdgeId>> Propagate(EdgeAnnotation& edge_annotation) const {
        std::map<bin_id, set<EdgeId>> answer;
        DEBUG("Propagating with propagator: " << name_);
        for (bin_id bin : edge_annotation.interesting_bins()) {
            DEBUG("Processing bin " << bin << " with propagator: " << name_);
            auto init_edges = edge_annotation.EdgesOfBin(bin, edge_length_threshold_);
            DEBUG("Initial edge cnt " << init_edges.size() << " (edge length threshold " << edge_length_threshold_ << ")");
            auto raw_propagated = PropagateEdges(init_edges);
            set<EdgeId> propagated;
            std::set_difference(raw_propagated.begin(), raw_propagated.end(),
                                init_edges.begin(), init_edges.end(),
                                std::inserter(propagated, propagated.end()));
            answer[bin] = std::move(propagated);
        }
        DEBUG("Finished propagating with propagator: " << name_);
        return answer;
    }

    virtual ~EdgeAnnotationPropagator() {}
private:
    DECL_LOGGER("EdgeAnnotationPropagator");
};

class ConnectingPathPropagator : public EdgeAnnotationPropagator {
    size_t path_length_threshold_;
    size_t path_edge_cnt_;
    const EdgeAnnotation& debug_annotation_;

    bin_id DetermineBin(const set<EdgeId>& edges) const {
        map<bin_id, size_t> cnt_map;
        for (EdgeId e : edges) {
            for (auto b : debug_annotation_.Annotation(e)) {
                cnt_map[b]++;
            }
        }
        bin_id candidate = "";
        for (auto cnt_el : cnt_map) {
            if (cnt_el.second > edges.size() / 2) {
                if (candidate.empty())
                    candidate = cnt_el.first;
                else
                    return "";
            }
        }
        return candidate;
    }

    bool BadPath(const vector<EdgeId>& path, bin_id base_bin) const {
        size_t cnt = 0;
        for (EdgeId e : path) {
            if (g().length(e) < 2000) 
                continue;
            auto ann = debug_annotation_.Annotation(e);
            if (!ann.empty() &&
                std::find(ann.begin(), ann.end(), base_bin) == ann.end()) {
                cnt++;
            }
        }
        return cnt > 0;
    }

    set<VertexId> CollectEdgeStarts(const set<EdgeId>& edges) const {
        set<VertexId> answer;
        for (EdgeId e : edges) {
            answer.insert(g().EdgeStart(e));
        }
        return answer;
    }

    set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const override {
        //static size_t pic_cnt = 0;
        bin_id bin = DetermineBin(edges);
        if (!bin.empty()) {
            DEBUG("Bin determined as " << bin);
        } else {
            DEBUG("Failed to determine bin");
        }
        set<EdgeId> answer;
        set<VertexId> starts = CollectEdgeStarts(edges);
        for (EdgeId e : edges) {
            PathProcessor<Graph> path_searcher(g(), g().EdgeEnd(e), path_length_threshold_);
            for (VertexId v : starts) {
                auto callback = AdapterCallback<Graph>([&](const vector<EdgeId>& path) {
                    //if (pic_cnt < 10) {
                    //if (BadPath(path, bin)) {
                    //    auto to_draw = path;
                    //    to_draw.insert(to_draw.begin(), e);
                    //    PrintAnnotatedAlongPath(gp(), to_draw, debug_annotation_, "/home/snurk/tmp/pics/pic_" + ToString(++pic_cnt) + "_");
                    //}
                    //}
                    insert_all(answer, path);
                }, true);
                TRACE("Launching path search between edge " << g().str(e) << " and vertex "
                        << g().str(v) << " with length bound " << path_length_threshold_);
                path_searcher.Process(v, 0, path_length_threshold_, callback, path_edge_cnt_);
            }
        }
        return answer;
    }

public:
    ConnectingPathPropagator(const conj_graph_pack& gp,
                             size_t path_length_threshold, 
                             size_t path_edge_cnt,
                             const EdgeAnnotation& ann) :
        EdgeAnnotationPropagator(gp, "ConnectingPath"),
        path_length_threshold_(path_length_threshold),
        path_edge_cnt_(path_edge_cnt),
        debug_annotation_(ann) {}

private:
    DECL_LOGGER("ConnectingPathPropagator");
};

//FIXME make threshold coverage-aware
class PairedInfoPropagator : public EdgeAnnotationPropagator {
    omnigraph::de::DEWeight weight_threshold_;
    set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const override {
        set<EdgeId> answer;
        for (EdgeId e1 : edges) {
            DEBUG("Searching for paired neighbours of " << g().str(e1));
            for (const auto& index : gp().clustered_indices)
                for (auto i : index.Get(e1))
                    for (auto point : i.second)
                        if (math::ge(point.weight, weight_threshold_)) {
                            DEBUG("Adding (" << g().str(e1) << "," << g().str(i.first) << "); " << point);
                            answer.insert(i.first);
                        }	    
        }
        return answer;
    }
public:
    PairedInfoPropagator(const conj_graph_pack& gp, omnigraph::de::DEWeight threshold):
        EdgeAnnotationPropagator(gp, "PairedInfo"), weight_threshold_(threshold) {}
private:
    DECL_LOGGER("PairedInfoPropagator");
};

class ContigPropagator : public EdgeAnnotationPropagator {
public:
    ContigPropagator(const conj_graph_pack& gp,
                     io::SingleStream& contigs) :
        EdgeAnnotationPropagator(gp, "ContigPropagator"),
        contigs_(contigs),
        mapper_(MapperInstance(gp))
    {}
protected:
    set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const override {
        contigs_.reset();
        set<EdgeId> answer;
        io::SingleRead contig;
        while (!contigs_.eof()) {
            contigs_ >> contig;
            auto edges_of_contig = mapper_->MapRead(contig).simple_path();
            for (EdgeId e : edges_of_contig) {
                if (edges.count(e)) {
                    DEBUG("Edge " << gp().g.str(e) << " belongs to the contig #" << 
                            contig.name() << " of " << edges_of_contig.size() << " edges");
                    insert_all(answer, edges_of_contig);
                    break;
                }
            }
        }
        return answer;
    }

private:
    io::SingleStream& contigs_;
    shared_ptr<SequenceMapper<Graph>> mapper_;

    DECL_LOGGER("ContigPropagator");
};

class TipPropagator : public EdgeAnnotationPropagator {

public:
    TipPropagator(const conj_graph_pack& gp) :
        EdgeAnnotationPropagator(gp, "TipPropagator"), tipper_(gp.g) {}

protected:
    set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const override {
        set<EdgeId> answer;
        for (EdgeId e1 : edges) {
            auto v = g().EdgeEnd(e1);
            auto neighbours = g().OutgoingEdges(v);
            auto e2_it = std::find_if(neighbours.begin(), neighbours.end(), [&](EdgeId e2){return edges.count(e2);});
            if (e2_it == neighbours.end()) {
                TRACE(e1.int_id() << " has no neighbours from the same bin");
                continue;
            }
            TRACE("Finding tips between " << e1.int_id() << " and " << e2_it->int_id());
            for (EdgeId posTip : g().IncidentEdges(v)) {
                if (edges.count(posTip))
                    continue;
                TRACE("Checking " << posTip.int_id() << "...");
                if (tipper_.Check(posTip)) {
                    TRACE("A tip is found!");
                    answer.insert(posTip);
                }
            }
        }
        return answer;
    }

private:
    TipCondition<Graph> tipper_;
    DECL_LOGGER("TipPropagator");
};

class AnnotationChecker {
    const Graph& g_;
    const EdgeAnnotation& edge_annotation_;
    size_t edge_length_threshold_;
public:
    AnnotationChecker(const Graph& g,
                      const EdgeAnnotation& edge_annotation,
                      size_t edge_length_threshold = EDGE_LENGTH_THRESHOLD) :
        g_(g),
        edge_annotation_(edge_annotation),
        edge_length_threshold_(edge_length_threshold) {
    }

    size_t Check(bin_id bin, const set<EdgeId>& propagated_edges) {
        DEBUG("Checking edges to be annotated with " << bin);
        size_t answer = 0;
        for (EdgeId e : propagated_edges) {
            if (g_.length(e) < edge_length_threshold_)
                continue;
            auto ann = edge_annotation_.Annotation(e);
            for (auto b : ann) {
                if (b != bin) {
                    DEBUG("Edge " << g_.str(e) << " already was annotated as " << b);
                    ++answer;
                    break;
                }
            }
        }
        return answer;
    }

private:
    DECL_LOGGER("AnnotationChecker");
};

void AnnotationPropagator::Run(io::SingleStream& /*contigs*/, 
                     EdgeAnnotation& edge_annotation
                     /*const string& annotation_out_fn*/) {
    std::vector<std::shared_ptr<EdgeAnnotationPropagator>> propagator_pipeline {
        std::make_shared<ConnectingPathPropagator>(gp_, 8000, 10, edge_annotation),
        std::make_shared<TipPropagator>(gp_), 
        std::make_shared<PairedInfoPropagator>(gp_, 10.)};//,
//        std::make_shared<ContigPropagator>(gp_, contigs)};//,
//        std::make_shared<ConnectingPathPropagator>(gp_, 8000, 10, edge_annotation),
//        std::make_shared<ContigPropagator>(gp_, contigs),
//        std::make_shared<TipPropagator>(gp_)};

    AnnotationChecker checker(gp_.g, edge_annotation);

    for (const auto& bin_id : edge_annotation.interesting_bins()) {
        size_t problem_cnt = checker.Check(bin_id, edge_annotation.EdgesOfBin(bin_id, EDGE_LENGTH_THRESHOLD));
        DEBUG("Bin " << bin_id << " had " << problem_cnt << " problems");
    }

    for (auto prop_ptr : propagator_pipeline) {
        DEBUG("Propagating with: " << prop_ptr->name());
        auto propagation_map = prop_ptr->Propagate(edge_annotation);

        DEBUG("Extending " << propagation_map.size() << " bins after propagation with: " << prop_ptr->name());
        for (const auto& bin_prop : propagation_map) {
            const auto& bin_id = bin_prop.first;
            const auto& edges = bin_prop.second;
            DEBUG("Extending bin " << bin_id << " with "
                      << edges.size() << " edges and their conjugates");
            size_t problem_cnt = checker.Check(bin_id, edges);
            DEBUG("Propagation of bin " << bin_id << " with " << prop_ptr->name()
                      << " lead to " << problem_cnt << " binning problems");
            edge_annotation.StickAnnotation(edges, bin_id);
        }
        DEBUG("Applied bin extensions from propagator " << prop_ptr->name());
    }
}

}
