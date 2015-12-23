//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "annotation.hpp"
#include "simple_tools.hpp"
#include "logger/log_writers.hpp"

#include "graphio.hpp"
#include "graph_pack.hpp"
#include "io/file_reader.hpp"

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

    EdgeAnnotation Propagate(const EdgeAnnotation& edge_annotation) const {
        EdgeAnnotation extended_annotation(edge_annotation);
        DEBUG("Propagating");
        for (bin_id bin : edge_annotation.interesting_bins()) {
            DEBUG("Processing bin " << bin);
            extended_annotation.StickAnnotation(PropagateEdges(edge_annotation.EdgesOfBin(bin)), bin);
        }
        return extended_annotation;
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
        set<EdgeId> answer(edges);
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

//todo add unbranching reach and pair-info aware propagation

class AnnotationPropagator {
    const conj_graph_pack& gp_;

    void DumpContigAnnotation(io::SingleStream& contigs,
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

public:
    AnnotationPropagator(const conj_graph_pack& gp) :
                     gp_(gp) {
    }

    void Run(io::SingleStream& contigs, const string& annotation_in_fn,
                         const vector<bin_id>& bins_of_interest, 
                         const string& annotation_out_fn) {
        AnnotationStream annotation_in(annotation_in_fn);
        EdgeAnnotation edge_annotation(gp_, bins_of_interest);
        edge_annotation.Fill(contigs, annotation_in);
        //FIXME magic constants
        size_t iteration_cnt = 2;
        ConnectingPathPropagator edge_propagator(gp_, 6000);

        auto annotation_ptr = make_shared<EdgeAnnotation>(edge_annotation);
        for (size_t i = 0; i < iteration_cnt; ++i) {
             annotation_ptr = make_shared<EdgeAnnotation>(
                                             edge_propagator.Propagate(*annotation_ptr));
        }

        contigs.reset();
        DumpContigAnnotation(contigs, *annotation_ptr, annotation_out_fn);
    }
};

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
    string saves_path = argv[2];
    string contigs_path = argv[3];
    string annotation_in_fn = argv[4];
    string annotation_out_fn = argv[5];
//    debruijn_graph::Launch(k, saves_path, contigs_path);

    std::vector<bin_id> bins_of_interest;
    for (int i = 6; i < argc; ++i) {
        bins_of_interest.push_back(argv[i]);
    }

    conj_graph_pack gp(k, "tmp", 0);
    gp.kmer_mapper.Attach();
    INFO("Load graph from " << saves_path);
    graphio::ScanGraphPack(saves_path, gp);
    auto contigs_stream_ptr = make_shared<io::FileReadStream>(contigs_path);

    AnnotationPropagator propagator(gp);
    propagator.Run(*contigs_stream_ptr, annotation_in_fn, bins_of_interest, annotation_out_fn);

    return 0;
}
