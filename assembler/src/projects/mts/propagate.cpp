//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dev_support/simple_tools.hpp"
#include "dev_support/logger/log_writers.hpp"

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
#include "io/reads_io/file_reader.hpp"
#include "algorithms/simplification/tip_clipper.hpp"
#include "propagate.hpp"

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

namespace debruijn_graph {

class EdgeAnnotationPropagator {
    const conj_graph_pack& gp_;

protected:
    const conj_graph_pack& gp() const {
        return gp_;
    }

    const Graph& g() const {
        return gp_.g;
    }

    virtual set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const = 0;

public:
    EdgeAnnotationPropagator(const conj_graph_pack& gp) : gp_(gp) {}

    void Propagate(EdgeAnnotation& edge_annotation) const {
        DEBUG("Propagating");
        for (bin_id bin : edge_annotation.interesting_bins()) {
            auto edges = edge_annotation.EdgesOfBin(bin);
            insert_all(edges, PropagateEdges(edges));
            DEBUG("Processing bin " << bin);
            edge_annotation.StickAnnotation(edges, bin);
        }
    }

    virtual ~EdgeAnnotationPropagator() {}
private:
    DECL_LOGGER("EdgeAnnotationPropagator");
};

class ConnectingPathPropagator : public EdgeAnnotationPropagator {
    size_t path_length_threshold_;

    set<VertexId> CollectEdgeStarts(const set<EdgeId>& edges) const {
        set<VertexId> answer;
        for (EdgeId e : edges) {
            answer.insert(g().EdgeStart(e));
        }
        return answer;
    }

    set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const override {
        set<EdgeId> answer;
        set<VertexId> starts = CollectEdgeStarts(edges);
        for (EdgeId e : edges) {
            PathProcessor<Graph> path_searcher(g(), g().EdgeEnd(e), path_length_threshold_);
            for (VertexId v : starts) {
                auto callback = AdaptorCallback<Graph>([&](vector<EdgeId> path) {
                    insert_all(answer, path);
                });
                DEBUG("Launching path search between edge " << g().str(e) << " and vertex "
                        << g().str(v) << " with length bound " << path_length_threshold_);
                path_searcher.Process(v, 0, path_length_threshold_, callback);
            }
        }
        return answer;
    }

public:
    ConnectingPathPropagator(const conj_graph_pack& gp, size_t path_length_threshold) :
        EdgeAnnotationPropagator(gp), path_length_threshold_(path_length_threshold) {}

private:
    DECL_LOGGER("ConnectingPathPropagator");
};

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
        EdgeAnnotationPropagator(gp), weight_threshold_(threshold) {}
private:
    DECL_LOGGER("PairedInfoPropagator");
};

class ContigPropagator : public EdgeAnnotationPropagator {

public:
    ContigPropagator(const conj_graph_pack& gp,
                     io::SingleStream& contigs,
                     EdgeAnnotation& annotation) :
        EdgeAnnotationPropagator(gp),
        contigs_(contigs),
        annotation_(annotation)
    {}
protected:
    set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const override {
        contigs_.reset();
        set<EdgeId> answer;
        io::SingleRead contig;
        while (!contigs_.eof()) {
            contigs_ >> contig;
            auto edges_of_contig = annotation_.EdgesOfContig(contig);
            for (EdgeId e : edges_of_contig) {
                if (edges.count(e)) {
                    DEBUG(e << " belongs to the contig #" << contig.name());
                    insert_all(answer, edges_of_contig);
                    break;
                }
            }
        }
        return answer;
    }

private:
    io::SingleStream& contigs_;
    const EdgeAnnotation& annotation_;
    DECL_LOGGER("ContigPropagator");
};

class TipPropagator : public EdgeAnnotationPropagator {

public:
    TipPropagator(const conj_graph_pack& gp) :
        EdgeAnnotationPropagator(gp), tipper_(gp.g) {}

protected:
    set<EdgeId> PropagateEdges(const set<EdgeId>& edges) const override {
        set<EdgeId> answer;
        for (EdgeId e1 : edges) {
            auto v = g().EdgeEnd(e1);
            auto neighbours = g().OutgoingEdges(v);
            auto e2_it = std::find_if(neighbours.begin(), neighbours.end(), [&](EdgeId e2){return edges.count(e2);});
            if (e2_it == neighbours.end()) {
                DEBUG(e1.int_id() << " has no neighbours from the same bin");
                continue;
            }
            DEBUG("Finding tips between " << e1.int_id() << " and " << e2_it->int_id());
            for (EdgeId posTip : g().IncidentEdges(v)) {
                if (edges.count(posTip))
                    continue;
                DEBUG("Checking " << posTip.int_id() << "...");
                if (tipper_.Check(posTip)) {
                    DEBUG("A tip is found!");
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
    ConnectingPathPropagator edge_propagator(gp_, 6000);
    TipPropagator tip_propagator(gp_);
    ContigPropagator contig_propagator(gp_, contigs, edge_annotation);
    PairedInfoPropagator paired_propagator(gp_, 2.0);

    //TODO: make this configurable
    std::vector<EdgeAnnotationPropagator*> propagator_pipeline =
        {&edge_propagator, &tip_propagator, &paired_propagator,
         &edge_propagator, &contig_propagator, &tip_propagator};

    for (auto prop_ptr : propagator_pipeline) {
        prop_ptr->Propagate(edge_annotation);
    }

    contigs.reset();
    DumpContigAnnotation(contigs, edge_annotation, annotation_out_fn);
}

}

int main(int argc, char** argv) {
    using namespace debruijn_graph;

    if (argc < 6) {
        cout << "Usage: annotation_propagator <K> <saves path> <contigs_path> "
                "<init binning info> <final binning info> (<bins of interest>)*"  << endl;
        exit(1);
    }
    //TmpFolderFixture fixture("tmp");
    create_console_logger();
    size_t k = lexical_cast<size_t>(argv[1]);
    string config_path = argv[2];
    string saves_path = argv[3];
    string contigs_path = argv[4];
    string annotation_in_fn = argv[5];
    string annotation_out_fn = argv[6];
//    debruijn_graph::Launch(k, saves_path, contigs_path);

    std::vector<bin_id> bins_of_interest;
    for (int i = 7; i < argc; ++i) {
        bins_of_interest.push_back(argv[i]);
    }

    cfg::create_instance(config_path);
    conj_graph_pack gp(k, "tmp", cfg::get().ds.reads.lib_count());
    gp.kmer_mapper.Attach();
    INFO("Load graph from " << saves_path);
    graphio::ScanWithClusteredIndices(saves_path, gp, gp.clustered_indices);
    auto contigs_stream_ptr = make_shared<io::FileReadStream>(contigs_path);

    AnnotationPropagator propagator(gp);
    propagator.Run(*contigs_stream_ptr, annotation_in_fn, bins_of_interest, annotation_out_fn);

    return 0;
}
