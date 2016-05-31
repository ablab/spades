//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dev_support/simple_tools.hpp"

//#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
//#include "io/reads_io/file_reader.hpp"
#include "algorithms/simplification/tip_clipper.hpp"
#include "propagate.hpp"
#include "visualization.hpp"

namespace debruijn_graph {

//FIXME 2kb edge length threshold might affect tip clipper in undesired way
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

//    size_t CountWrongCAG(EdgeAnnotation& edge_annotation, const set<EdgeId>& edges, bin_id bin) const {
//        bin_id to_check = (bin == "CAG1") ? "CAG2" : "CAG1";
//        size_t answer = 0;
//        for (EdgeId e : edges) {
//            auto ann = edge_annotation.Annotation(e);
//            if (std::find(ann.begin(), ann.end(), to_check) != ann.end()) {
//                if (g().length(e) > edge_length_threshold_)
//                    ++answer;
//            }
//        }
//        return answer;
//    }

public:
    EdgeAnnotationPropagator(const conj_graph_pack& gp,
                             const string& name,
                             size_t edge_length_threshold = 2000) :
                    gp_(gp),
                    name_(name),
                    edge_length_threshold_(edge_length_threshold) {}

    void Propagate(EdgeAnnotation& edge_annotation) const {
        DEBUG("Propagating with propagator: " << name_);
        for (bin_id bin : edge_annotation.interesting_bins()) {
            DEBUG("Processing bin " << bin << " with propagator: " << name_);
            auto edges = edge_annotation.EdgesOfBin(bin, edge_length_threshold_);
            size_t init_cnt = edges.size();
            DEBUG("Initial propagation edge cnt " << init_cnt << " (edge length threshold " << edge_length_threshold_ << ")");
//            DEBUG("WRONG CNT " << CountWrongCAG(edge_annotation, edges, bin));
            insert_all(edges, PropagateEdges(edges));
            DEBUG("Propagated on " << (edges.size() - init_cnt) << " edges");
//            DEBUG("WRONG CNT " << CountWrongCAG(edge_annotation, edges, bin));
            DEBUG("Sticking annotation to edges and conjugates");
            edge_annotation.StickAnnotation(edges, bin);
            DEBUG("Post-propagation bin edge cnt " << edge_annotation.EdgesOfBin(bin).size());
//            DEBUG("WRONG CNT " << CountWrongCAG(edge_annotation, edge_annotation.EdgesOfBin(bin), bin));
        }
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
        static size_t pic_cnt = 0;
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
            TRACE("Searching for paired neighbours of " << g().str(e1));
            for (const auto& index : gp().clustered_indices)
                for (auto i : index.Get(e1))
                    for (auto point : i.second)
                        if (math::ge(point.weight, weight_threshold_)) {
                            TRACE("Adding (" << g().str(e1) << "," << g().str(i.first) << "); " << point);
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
                    TRACE(e << " belongs to the contig #" << contig.name());
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

void AnnotationPropagator::DumpContigAnnotation(io::SingleStream& contigs,
                          const EdgeAnnotation& annotation,
                          const string& annotation_out_fn) const {
    AnnotationOutStream annotation_out(annotation_out_fn);
    io::SingleRead contig;
    while (!contigs.eof()) {
        contigs >> contig;
        auto relevant_bins = annotation.RelevantBins(contig);
        if (!relevant_bins.empty()) {
            annotation_out << ContigAnnotation(GetId(contig),
                                vector<bin_id>(relevant_bins.begin(), relevant_bins.end()));
        }
    }
}

void AnnotationPropagator::Run(io::SingleStream& contigs, const string& annotation_in_fn,
                     const vector<bin_id>& bins_of_interest,
                     const string& annotation_out_fn) {
    AnnotationStream annotation_in(annotation_in_fn);
    EdgeAnnotation edge_annotation(gp_, bins_of_interest);
    edge_annotation.Fill(contigs, annotation_in);

    //TODO: make this configurable
    std::vector<std::shared_ptr<EdgeAnnotationPropagator>> propagator_pipeline {
        std::make_shared<ConnectingPathPropagator>(gp_, 8000, 10, edge_annotation),
        std::make_shared<TipPropagator>(gp_), 
        std::make_shared<PairedInfoPropagator>(gp_, 10.),
        std::make_shared<ContigPropagator>(gp_, contigs)};//,
//        std::make_shared<ConnectingPathPropagator>(gp_, 8000, 10, edge_annotation),
//        std::make_shared<ContigPropagator>(gp_, contigs),
//        std::make_shared<TipPropagator>(gp_)};

    for (auto prop_ptr : propagator_pipeline) {
        prop_ptr->Propagate(edge_annotation);
    }

    contigs.reset();
    DumpContigAnnotation(contigs, edge_annotation, annotation_out_fn);
}

}
