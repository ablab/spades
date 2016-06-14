/*
 * stats.cpp
 *
 *  Created on: 3 Dec 2015
 *      Author: idmit
 */

#include "pipeline/graphio.hpp"
#include "pipeline/graph_pack.hpp"
#include "dev_support/simple_tools.hpp"
#include "dev_support/path_helper.hpp"
#include "math/xmath.h"
#include <iostream>
#include <vector>
#include "io/reads_io/multifile_reader.hpp"
#include "io/reads_io/splitting_wrapper.hpp"
#include "io/reads_io/modifying_reader_wrapper.hpp"
#include "io/reads_io/vector_reader.hpp"
#include "io/reads_io/file_reader.hpp"
#include "annotation.hpp"
#include "visualization.hpp"
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

    AnnotationStream annotation_stream(annotation_path);
    edge_annotation.Fill(contigs_stream, annotation_stream);
    return edge_annotation;
}

int main(int argc, char** argv) {
    using namespace GetOpt;

    size_t k;
    string saves_path, contigs_path, edges_path;
    vector<string> genomes_path;
    string annotation_in_fn, prop_annotation_in_fn;
    string graph_dir;
    vector<bin_id> bins_of_interest;

    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('s', saves_path)
            >> Option('r', genomes_path)
            >> Option('c', contigs_path)
            >> Option('a', annotation_in_fn)
            >> Option('e', edges_path)
            >> Option('p', prop_annotation_in_fn)
            //>> Option('o', graph_dir, "")
            >> Option('b', bins_of_interest, {})
        ;
    } catch(GetOptEx &ex) {
        cout << "Usage: stats -k <K> -s <saves path> -r <genomes path>+ "
                "-c <contigs_path> -a <init binning info> -e <edges_path> -p <propagated binning info> "
                "[-o <graph directory> (currently disabled)] [-b (<bins of interest>)+]"
             << endl;
        exit(1);
    }
    //TmpFolderFixture fixture("tmp");

    conj_graph_pack gp(k, "tmp", 0);
    gp.kmer_mapper.Attach();
    INFO("Load graph from " << saves_path);
    graphio::ScanGraphPack(saves_path, gp);
    gp.edge_pos.Attach();

    std::cout << "Reference\tEdges in genome alignment\t"
              << "Length in genome alignment\t"
              << "Unbinned edges\t"
              << "Length of unbinned edges\t"
              << "Propagated edges\t"
              << "Length of propagated edges\t"
              << endl;

    for (const auto genome_path : genomes_path) {
        auto ref_name = path::basename(genome_path);
        io::SingleRead genome = ReadGenome(genome_path);

        FillPos(gp, genome_path, "", true);

        io::FileReadStream contigs_stream(contigs_path);
        EdgeAnnotation edge_annotation = LoadAnnotation(
            gp, bins_of_interest, contigs_stream, annotation_in_fn);

        io::FileReadStream edges_stream(edges_path);
        EdgeAnnotation prop_edge_annotation = LoadAnnotation(
            gp, bins_of_interest, edges_stream, prop_annotation_in_fn);

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
                    /*std::cout << e.int_id() << "\t"
                              << gp.g.length(e) << "\t"
                              << range.size() << std::endl;*/
                    if (!graph_dir.empty()) {
                        std::string dot_export_path =
                            graph_dir + "/" + ref_name + "/" + std::to_string(e.int_id()) + ".dot";
                        PrintColoredAnnotatedGraphAroundEdge(
                            gp, e, prop_edge_annotation, dot_export_path);
                    }
                }
            }
        }

        cout << ref_name << "\t"
             << total_edge_cntr << "\t"
             << total_edgelen_cntr << "\t"
             << unbinned_edge_cntr << "\t"
             << unbinned_edgelen_cntr << "\t"
             << prop_binned_edge_cntr << "\t"
             << prop_binned_edgelen_cntr << "\t"
             << endl;
    }
}
