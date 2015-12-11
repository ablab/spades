/*
 * stats.cpp
 *
 *  Created on: 3 Dec 2015
 *      Author: idmit
 */

#include "graphio.hpp"
#include "graph_pack.hpp"
#include "simple_tools.hpp"
#include "debruijn_graph.hpp"
#include "xmath.h"
#include <iostream>
#include <vector>
#include "io/multifile_reader.hpp"
#include "io/splitting_wrapper.hpp"
#include "io/modifying_reader_wrapper.hpp"
#include "io/vector_reader.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "io/file_reader.hpp"
#include "annotation.hpp"
#include "omni/visualization/graph_colorer.hpp"

using namespace debruijn_graph;

io::SingleRead ReadSequence(io::SingleStream& reader) {
    VERIFY(!reader.eof());
    io::SingleRead read;
    reader >> read;
    return read;
}

io::SingleRead ReadGenome(const string& genome_path) {
    path::CheckFileExistenceFATAL(genome_path);
    auto genome_stream_ptr = std::make_shared<io::FileReadStream>(genome_path);
    return ReadSequence(*genome_stream_ptr);
}

EdgeAnnotation LoadAnnotation(const conj_graph_pack& gp,
                              const vector<bin_id>& bins_of_interest,
                              io::SingleStream* contigs_stream_ptr,
                              const string& annotation_path) {
    EdgeAnnotation edge_annotation(gp, bins_of_interest);
    EdgeAnnotation prop_edge_annotation(gp, bins_of_interest);

    AnnotationStream annotation_stream(annotation_path);
    edge_annotation.Fill(*contigs_stream_ptr, annotation_stream);
    return edge_annotation;
}

template <class Graph>
class AnnotatedGraphColorer
    : public omnigraph::visualization::GraphColorer<Graph> {
   public:
    AnnotatedGraphColorer(const EdgeAnnotation& annotation)
        : _annotation(annotation) {}

    string GetValue(typename Graph::VertexId) const { return "black"; }

    string GetValue(typename Graph::EdgeId edge) const {
        if (_annotation.Annotation(edge).empty()) {
            return "black";
        }
        return "red";
    }

   private:
    EdgeAnnotation _annotation;
};

void PrintColoredAnnotatedGraphAroundEdge(const conj_graph_pack& gp,
                                          const EdgeId& edge,
                                          const EdgeAnnotation& annotation,
                                          const string& output_filename) {
    std::cout << output_filename << std::endl;
    DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    auto colorer_ptr =
        std::make_shared<AnnotatedGraphColorer<Graph>>(annotation);
    GraphComponent<Graph> component = omnigraph::EdgeNeighborhood(gp.g, edge);
    omnigraph::visualization::WriteComponent<Graph>(component, output_filename,
                                                    colorer_ptr, labeler);
}

int main(int argc, char** argv) {
    if (argc < 7) {
        cout << "Usage: annotation_propagator <K> <saves path> <genome path> "
                "<contigs_path>  <init binning info> <propagated binning info> "
                "(<bins of interest>)*"
             << endl;
        exit(1);
    }
    TmpFolderFixture fixture("tmp");
    //    create_console_logger();
    size_t k = lexical_cast<size_t>(argv[1]);
    string saves_path = argv[2];
    string genome_path = argv[3];
    string contigs_path = argv[4];
    string annotation_in_fn = argv[5];
    string prop_annotation_in_fn = argv[6];

    std::vector<bin_id> bins_of_interest;
    for (int i = 6; i < argc; ++i) {
        bins_of_interest.push_back(argv[i]);
    }

    shared_ptr<io::SingleStream> contigs_stream_ptr =
        make_shared<io::FileReadStream>(contigs_path);

    io::SingleRead genome = ReadGenome(genome_path);

    conj_graph_pack gp(k, "tmp", 0);
    gp.kmer_mapper.Attach();
    INFO("Load graph from " << saves_path);
    graphio::ScanGraphPack(saves_path, gp);

    EdgeAnnotation edge_annotation = LoadAnnotation(
        gp, bins_of_interest, contigs_stream_ptr.get(), annotation_in_fn);
    contigs_stream_ptr->reset();
    EdgeAnnotation prop_edge_annotation = LoadAnnotation(
        gp, bins_of_interest, contigs_stream_ptr.get(), prop_annotation_in_fn);

    shared_ptr<SequenceMapper<Graph>> mapper(MapperInstance(gp));

    std::size_t prop_binned_edge_cntr = 0;
    std::size_t prop_binned_edgelen_cntr = 0;

    std::size_t unbinned_edge_cntr = 0;
    std::size_t unbinned_edgelen_cntr = 0;

    std::size_t total_edge_cntr = 0;
    std::size_t total_edgelen_cntr = 0;

    auto genome_graph_path = mapper->MapRead(genome).simple_path();

    gp.EnsurePos();
    for (EdgeId e : genome_graph_path) {
        total_edge_cntr++;
        total_edgelen_cntr += gp.g.length(e);
        if (edge_annotation.Annotation(e).empty() &&
            !prop_edge_annotation.Annotation(e).empty()) {
            prop_binned_edge_cntr++;
            prop_binned_edgelen_cntr += gp.g.length(e);
        } else if (edge_annotation.Annotation(e).empty() &&
                   prop_edge_annotation.Annotation(e).empty()) {
            // Only check for prop_annotation is necessary
            unbinned_edge_cntr++;
            unbinned_edgelen_cntr += gp.g.length(e);
            std::cout << "Edge id: " + std::to_string(e.int_id()) << std::endl;
            std::string dot_export_path("/Users/idmit/edge");
            dot_export_path += std::to_string(e.int_id()) + ".dot";
            PrintColoredAnnotatedGraphAroundEdge(gp, e, prop_edge_annotation,
                                                 dot_export_path);
        }
    }

    std::cout << "Total number of edges in genome alignment: "
              << total_edge_cntr << std::endl;
    std::cout << "Total length of edges in genome alignment: "
              << total_edgelen_cntr << std::endl;

    std::cout << "Total number of unbinned edges: " << unbinned_edge_cntr
              << std::endl;
    std::cout << "Total length of unbinned edges: " << unbinned_edgelen_cntr
              << std::endl;

    std::cout << "Total number of edges binned due to propagation: "
              << prop_binned_edge_cntr << std::endl;
    std::cout << "Total length of edges binned due to propagation: "
              << prop_binned_edgelen_cntr << std::endl;
}
