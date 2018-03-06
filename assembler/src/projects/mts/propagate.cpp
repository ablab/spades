//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/stl_utils.hpp"

//#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
//#include "io/reads_io/file_reader.hpp"
#include "modules/simplification/tip_clipper.hpp"
#include "propagate.hpp"
#include "visualization.hpp"

namespace debruijn_graph {
static const size_t EDGE_LENGTH_THRESHOLD = 2000,
                    EDGE_UPPER_THRESHOLD = 3000;

//FIXME 2kb edge length threshold might affect tip propagator in undesired way
class EdgeAnnotationPropagator {
    const conj_graph_pack& gp_;
    const string name_;
    size_t edge_length_threshold_, edge_upper_threshold_;

protected:
    size_t edge_length_threshold() const {
        return edge_length_threshold_;
    }

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
                             size_t edge_length_threshold = EDGE_LENGTH_THRESHOLD,
                             size_t edge_upper_threshold = EDGE_UPPER_THRESHOLD) :
                    gp_(gp),
                    name_(name),
                    edge_length_threshold_(edge_length_threshold),
                    edge_upper_threshold_(edge_upper_threshold) {}

    const std::string& name() const {
        return name_;
    }

    std::map<bin_id, set<EdgeId>> Propagate(EdgeAnnotation& edge_annotation) const {
        std::map<bin_id, set<EdgeId>> answer;
        DEBUG("Propagating with propagator: " << name_);
        for (bin_id bin : edge_annotation.interesting_bins()) {
            DEBUG("Processing bin " << bin << " with propagator: " << name_);
            auto init_edges = edge_annotation.EdgesOfBin(bin, edge_length_threshold());
            DEBUG("Initial edge cnt " << init_edges.size() << " (edge length threshold " << edge_length_threshold() << ")");
            auto raw_propagated = PropagateEdges(init_edges);
            auto old_size = raw_propagated.size();
            //Filter
            size_t n = 0;
            for (auto i = raw_propagated.begin(); i != raw_propagated.end(); ++n) {
                DEBUG("Edge cnt: " << raw_propagated.size() << "; iter " << n);
                if (gp_.g.length(*i) > edge_upper_threshold_)
                    raw_propagated.erase(i++);
                else
                    ++i;
            }
            DEBUG("Excluded " << (old_size - raw_propagated.size()) << " >" << edge_upper_threshold_ << "bp edges");
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

    set<VertexId> CollectEdgeStarts(const set<EdgeId>& edges) const {
        set<VertexId> answer;
        for (EdgeId e : edges) {
            answer.insert(g().EdgeStart(e));
        }
        return answer;
    }

    set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const override {
        DEBUG(__FUNCTION__);
        set<EdgeId> answer;
        set<VertexId> starts = CollectEdgeStarts(edges);
        for (EdgeId e : edges) {
            PathProcessor<Graph> path_searcher(g(), g().EdgeEnd(e), path_length_threshold_);
            for (VertexId v : starts) {
                auto callback = AdapterCallback<Graph>([&](const vector<EdgeId>& path) {
                    utils::insert_all(answer, path);
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
                             size_t length_threshold,
                             size_t path_length_threshold,
                             size_t path_edge_cnt) :
        EdgeAnnotationPropagator(gp, "ConnectingPath", length_threshold),
        path_length_threshold_(path_length_threshold),
        path_edge_cnt_(path_edge_cnt) {}

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
                for (auto i : index.Get(e1)) {
                    if (e1 == i.first) //No sense in self-propagation
                        continue;
                    for (auto point : i.second)
                        if (math::ge(point.weight, weight_threshold_)) {
                            TRACE("Adding (" << g().str(e1) << "," << g().str(i.first) << "); " << point);
                            answer.insert(i.first);
                        }
                }
        }
        return answer;
    }
public:
    PairedInfoPropagator(const conj_graph_pack& gp,
                         size_t length_threshold,
                         size_t upper_threshold,
                         omnigraph::de::DEWeight weight_threshold):
        EdgeAnnotationPropagator(gp, "PairedInfo", length_threshold, upper_threshold),
        weight_threshold_(weight_threshold)
    {}
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
                    utils::insert_all(answer, edges_of_contig);
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
    TipPropagator(const conj_graph_pack& gp, size_t length_threshold) :
        EdgeAnnotationPropagator(gp, "TipPropagator", length_threshold), tipper_(gp.g) {}

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
        MakePropagator<ConnectingPathPropagator>(8000, 10),
//        make_propagator<TipPropagator>(),
        MakePropagator<PairedInfoPropagator>(100500, 10.)};//,
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
