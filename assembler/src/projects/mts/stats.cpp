/*
 * stats.cpp
 *
 *  Created on: 3 Dec 2015
 *      Author: idmit
 */

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
#include "dev_support/simple_tools.hpp"
#include "math/xmath.h"
#include <iostream>
#include <vector>
#include "io/reads_io/multifile_reader.hpp"
#include "io/reads_io/splitting_wrapper.hpp"
#include "io/reads_io/modifying_reader_wrapper.hpp"
#include "io/reads_io/vector_reader.hpp"
//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/xml_parser.hpp>
#include "io/reads_io/file_reader.hpp"
#include "annotation.hpp"
#include "visualization/graph_colorer.hpp"
#include "visualization/position_filler.hpp"
#include "algorithms/simplification/tip_clipper.hpp"
#include "getopt_pp/getopt_pp.h"

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
                              io::SingleStream& contigs_stream,
                              const string& annotation_path) {
    EdgeAnnotation edge_annotation(gp, bins_of_interest);
    EdgeAnnotation prop_edge_annotation(gp, bins_of_interest);

    AnnotationStream annotation_stream(annotation_path);
    edge_annotation.Fill(contigs_stream, annotation_stream);
    return edge_annotation;
}

template <class Graph>
class AnnotatedGraphColorer
    : public omnigraph::visualization::GraphColorer<Graph> {
   public:
    AnnotatedGraphColorer(const EdgeAnnotation& annotation)
        : _annotation(annotation) {
        std::size_t i = 0;
        for (const auto& b_id : _annotation.interesting_bins()) {
            _color_map[b_id] = i;
            ++i;
        }
    }

    string GetValue(typename Graph::VertexId) const { return "black"; }

    string GetValue(typename Graph::EdgeId edge) const {
        if (_annotation.Annotation(edge).empty()) {
            return "black";
        }
        bin_id b_id = _annotation.Annotation(edge)[0];
        std::string color = _colors[_color_map[b_id]];
        for (std::size_t i = 1; i < _annotation.Annotation(edge).size(); ++i) {
            b_id = _annotation.Annotation(edge)[i];
            color += ":" + _colors[_color_map[b_id]];
        }
        return color;
    }

   private:
    EdgeAnnotation _annotation;
    static std::vector<std::string> _colors;
    mutable std::map<bin_id, std::size_t> _color_map;
};

template <class Graph>
std::vector<std::string> AnnotatedGraphColorer<Graph>::_colors = {
    "red", "blue", "green", "yellow", "orange", "purple", "pink"};

void PrintColoredAnnotatedGraphAroundEdge(const conj_graph_pack& gp,
                                          const EdgeId& edge,
                                          const EdgeAnnotation& annotation,
                                          const string& output_filename) {
    //std::cout << output_filename << std::endl;
    DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    auto colorer_ptr =
        std::make_shared<AnnotatedGraphColorer<Graph>>(annotation);
    GraphComponent<Graph> component = omnigraph::EdgeNeighborhood(gp.g, edge, 100, 10000);
    omnigraph::visualization::WriteComponent<Graph>(component, output_filename,
                                                    colorer_ptr, labeler);
}

int main(int argc, char** argv) {
    using namespace GetOpt;

    size_t k;
    string saves_path, genome_path, contigs_path;
    string annotation_in_fn, prop_annotation_in_fn;
    string graph_dir;
    std::vector<bin_id> bins_of_interest;

    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('s', saves_path)
            >> Option('r', genome_path)
            >> Option('c', contigs_path)
            >> Option('a', annotation_in_fn)
            >> Option('p', prop_annotation_in_fn)
            >> Option('o', graph_dir)
            >> Option('b', bins_of_interest, {})
        ;
    } catch(GetOptEx &ex) {
        cout << "Usage: stats -k <K> -s <saves path> -r <genome path> "
                "-c <contigs_path> -a <init binning info> -p <propagated binning info> "
                "-o <graph directory> [-b (<bins of interest>)+]"
             << endl;
        exit(1);
    }
    TmpFolderFixture fixture("tmp");

    io::FileReadStream contigs_stream(contigs_path);

    io::SingleRead genome = ReadGenome(genome_path);

    conj_graph_pack gp(k, "tmp", 0);
    gp.kmer_mapper.Attach();
    INFO("Load graph from " << saves_path);
    graphio::ScanGraphPack(saves_path, gp);
    gp.genome.SetSequence(genome.sequence());
    gp.edge_pos.Attach();
    FillPos(gp, gp.genome.GetSequence(), "ref0");
    FillPos(gp, !gp.genome.GetSequence(), "ref1");

    EdgeAnnotation edge_annotation = LoadAnnotation(
        gp, bins_of_interest, contigs_stream, annotation_in_fn);
    contigs_stream.reset();
    EdgeAnnotation prop_edge_annotation = LoadAnnotation(
        gp, bins_of_interest, contigs_stream, prop_annotation_in_fn);

    shared_ptr<SequenceMapper<Graph>> mapper(MapperInstance(gp));

    std::size_t prop_binned_edge_cntr = 0;
    std::size_t prop_binned_edgelen_cntr = 0;

    std::size_t unbinned_edge_cntr = 0;
    std::size_t unbinned_edgelen_cntr = 0;

    std::size_t total_edge_cntr = 0;
    std::size_t total_edgelen_cntr = 0;

    auto genome_graph_path = mapper->MapRead(genome);
    std::set<EdgeId> unbinned_edges;

    gp.EnsurePos();
    for (size_t i = 0; i < genome_graph_path.size(); ++i) {
        EdgeId e = genome_graph_path[i].first;
        auto range = genome_graph_path[i].second.mapped_range;
        total_edge_cntr++;
        total_edgelen_cntr += gp.g.length(e);
        if (edge_annotation.Annotation(e).empty() &&
            !prop_edge_annotation.Annotation(e).empty()) {
            DEBUG(e.int_id() << " was propagated\n");
            prop_binned_edge_cntr++;
            prop_binned_edgelen_cntr += gp.g.length(e);
        } else if (edge_annotation.Annotation(e).empty() &&
                   prop_edge_annotation.Annotation(e).empty()) {
            // Only check for prop_annotation is necessary
            if (unbinned_edges.count(e) == 0) {
                unbinned_edges.insert(e);
                unbinned_edge_cntr++;
                unbinned_edgelen_cntr += range.size();
                std::cout << e.int_id() << "\t"
                          << gp.g.length(e) << "\t"
                          << range.size() << std::endl;
                std::string dot_export_path(graph_dir);
                dot_export_path += std::to_string(e.int_id()) + ".dot";
                PrintColoredAnnotatedGraphAroundEdge(
                    gp, e, prop_edge_annotation, dot_export_path);
            }
        }
    }

    std::cerr << "Total number of edges in genome alignment: "
              << total_edge_cntr << std::endl;
    std::cerr << "Total length of edges in genome alignment: "
              << total_edgelen_cntr << std::endl;

    std::cerr << "Total number of unbinned edges: " << unbinned_edge_cntr
              << std::endl;
    std::cerr << "Total length of unbinned edges: " << unbinned_edgelen_cntr
              << std::endl;

    std::cerr << "Total number of edges binned due to propagation: "
              << prop_binned_edge_cntr << std::endl;
    std::cerr << "Total length of edges binned due to propagation: "
              << prop_binned_edgelen_cntr << std::endl;
}
