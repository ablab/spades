//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#include "annotation.hpp"
#include "visualization.hpp"

#include "getopt_pp/getopt_pp.h"
#include "io/binary/graph_pack.hpp"
#include "io/reads/file_reader.hpp"
#include "math/xmath.h"
#include "pipeline/graph_pack.hpp"
#include "pipeline/graph_pack_helpers.h"
#include "pipeline/sequence_mapper_gp_api.hpp"
#include "utils/stl_utils.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/logger/log_writers.hpp"
#include "visualization/position_filler.hpp"

#include <iostream>
#include <vector>

using namespace debruijn_graph;
using namespace std;

io::SingleRead ReadSequence(io::SingleStream& reader) {
    VERIFY(!reader.eof());
    io::SingleRead read;
    reader >> read;
    return read;
}

io::SingleRead ReadGenome(const filesystem::path& genome_path) {
    CHECK_FATAL_ERROR(exists(genome_path), "File " << genome_path << " doesn't exist or can't be read!");
    auto genome_stream_ptr = io::EasyStream(genome_path, false);
    return ReadSequence(genome_stream_ptr);
}

EdgeAnnotation LoadAnnotation(const graph_pack::GraphPack& gp,
                              const vector<bin_id>& bins_of_interest,
                              io::SingleStream& contigs_stream,
                              io::SingleStream& splits_stream,
                              const filesystem::path& annotation_path) {
    AnnotationFiller filler(gp, bins_of_interest);
    AnnotationStream annotation_stream(annotation_path);
    return filler(contigs_stream, splits_stream, annotation_stream);
}

class BinnedInfo : public pair<size_t, size_t> {
public:
    BinnedInfo(): pair(0, 0) {}
};

void add_edge_info(BinnedInfo& info, size_t edge_length) {
    ++info.first;
    info.second += edge_length;
}

ostream& operator<<(ostream& str, const BinnedInfo& info) {
    str << info.first << "\t" << info.second;
    return str;
}

void create_console_logger() {
    logging::logger *log = logging::create_logger("", logging::L_INFO);
    log->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(log);
}

int main(int argc, char** argv) {
    create_console_logger();

    using namespace GetOpt;

    size_t k;
    filesystem::path saves_path, contigs_path, splits_path, edges_path;
    vector<filesystem::path> genomes_path;
    filesystem::path annotation_in_fn, prop_annotation_in_fn;
    filesystem::path table_fn, graph_dir;
    vector<bin_id> bins_of_interest;

    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('s', saves_path)
            >> Option('r', genomes_path)
            >> Option('c', contigs_path)
            >> Option('f', splits_path)
            >> Option('a', annotation_in_fn)
            >> Option('e', edges_path)
            >> Option('p', prop_annotation_in_fn)
            >> Option('o', table_fn)
            //>> Option('d', graph_dir, "")
            >> Option('b', bins_of_interest, {})
        ;
    } catch(GetOptEx &ex) {
        cout << "Usage: stats -k <K> -s <saves path> -r <genomes path>+ "
                "-f <splits_path> -c <contigs_path> -a <init binning info> -e <edges_path> -p <propagated binning info> "
                "-o <stats table> [-d <graph directory> (currently disabled)] [-b (<bins of interest>)+]"
             << endl;
        exit(1);
    }
    //TmpFolderFixture fixture("tmp");

    graph_pack::GraphPack gp(k, "tmp", 0);
    gp.get_mutable<KmerMapper<Graph>>().Attach();
    INFO("Load graph from " << saves_path);
    io::binary::BasePackIO().Load(saves_path, gp);
    gp.get_mutable<EdgesPositionHandler<Graph>>().Attach();

    ofstream output(table_fn);

    output << "Reference\t"
           << "Aligned edges\tAlignment length\t"
           << "Binned edges\tBinned length\t"
           << "Unbinned edges\tUnbinned length\t"
           << "Pre-binned edges\tPre-binned length\t"
           << "Propagated edges\tPropagated length" << endl;

    for (const auto &genome_path : genomes_path) {
        auto ref_name = genome_path.stem();
        io::SingleRead genome = ReadGenome(genome_path);

        visualization::position_filler::FillPos(gp, genome_path, "", true);

        auto contigs_stream = io::EasyStream(contigs_path, false);
        auto splits_stream = io::EasyStream(splits_path, false);
        EdgeAnnotation edge_annotation = LoadAnnotation(
            gp, bins_of_interest, contigs_stream, 
            splits_stream, annotation_in_fn);

        auto edges_stream = io::EasyStream(edges_path, false);
        auto edges_stream2 = io::EasyStream(edges_path, false);
        EdgeAnnotation prop_edge_annotation = LoadAnnotation(
            gp, bins_of_interest, 
            edges_stream, edges_stream2, 
            prop_annotation_in_fn);

        shared_ptr<SequenceMapper<Graph>> mapper(MapperInstance(gp));

        BinnedInfo pre_binned_info, prop_binned_info, binned_info,
                   unbinned_info, total_info;

        auto genome_graph_path = mapper->MapRead(genome);
        std::set<EdgeId> unbinned_edges;

        auto const &graph = gp.get<Graph>();
        EnsurePos(gp);
        for (size_t i = 0; i < genome_graph_path.size(); ++i) {
            EdgeId e = genome_graph_path[i].first;
            auto range = genome_graph_path[i].second.mapped_range;
            add_edge_info(total_info, graph.length(e));
            if (edge_annotation.Annotation(e).empty()) {
                if (prop_edge_annotation.Annotation(e).empty()) {
                    // Only check for prop_annotation is necessary
                    if (unbinned_edges.count(e) == 0) {
                        unbinned_edges.insert(e);
                        add_edge_info(unbinned_info, range.size());
                        /*std::cout << e.int_id() << "\t"
                                  << gp.g.length(e) << "\t"
                                  << range.size() << std::endl;*/
                        if (!graph_dir.empty()) {
                            std::filesystem::path dot_export_path =
                                graph_dir / (ref_name / (std::to_string(e.int_id()) + ".dot"));
                            PrintColoredAnnotatedGraphAroundEdge(
                                gp, e, prop_edge_annotation, dot_export_path);
                        }
                    }
                } else {
                    DEBUG(e.int_id() << " was propagated\n");
                    add_edge_info(prop_binned_info, graph.length(e));
                    add_edge_info(binned_info, graph.length(e));
                }
            } else {
                add_edge_info(pre_binned_info, graph.length(e));
                if (prop_edge_annotation.Annotation(e).empty()) {
                    WARN(e.int_id() << " was lost during propagation\n");
                } else {
                    add_edge_info(binned_info, graph.length(e));
                }
            }
        }

        output << ref_name         << "\t"
               << total_info       << "\t"
               << binned_info      << "\t"
               << unbinned_info    << "\t"
               << pre_binned_info  << "\t"
               << prop_binned_info << endl;
    }
}
