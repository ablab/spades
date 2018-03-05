//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "getopt_pp/getopt_pp.h"
#include "io/reads/io_helper.hpp"
#include "io/reads/osequencestream.hpp"
#include "pipeline/graphio.hpp"
#include "logger.hpp"
#include "read_binning.hpp"
#include "propagate.hpp"
#include "visualization/position_filler.hpp"

using namespace debruijn_graph;

std::string add_suffix(const std::string& path, const std::string& suffix) {
    auto ext = fs::extension(path);
    return path.substr(0, path.length() - ext.length()) + suffix + ext;
}

//TODO: refactor to process the graph only once
void DumpEdges(const Graph& g, const string& out_edges) {
    INFO("Dumping edges to " << out_edges);
    io::OFastaReadStream oss(out_edges);
    for (auto it = g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        oss << io::SingleRead("NODE_" + std::to_string(g.int_id(e)), g.EdgeNucls(e).str());
    }
}

void DumpAnnotation(const Graph& g, const EdgeAnnotation& edge_annotation, const string& out_annotation) {
    INFO("Dumping annotation to " << out_annotation);
    AnnotationOutStream annotation_out(out_annotation);
    for (auto it = g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        auto relevant_bins = edge_annotation.Annotation(e);
        if (!relevant_bins.empty()) {
            annotation_out << ContigAnnotation("NODE_" + std::to_string(g.int_id(e)),
                                               vector<bin_id>(relevant_bins.begin(), relevant_bins.end()));
        }
    }
}

int main(int argc, char** argv) {
    using namespace GetOpt;

    //TmpFolderFixture fixture("tmp");
    create_console_logger();

    size_t k;
    string saves_path, contigs_path, splits_path, annotation_path, bins_file;
    vector<string> sample_names, left_reads, right_reads;
    string out_root, edges_dump, propagation_dump;
    size_t length_threshold;
    bool no_binning;
    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('s', saves_path)
            >> Option('c', contigs_path)
            >> Option('f', splits_path)
            >> Option('a', annotation_path)
            >> Option('t', length_threshold, (size_t)2000)
            >> Option('b', bins_file, "")
            >> Option('n', "names", sample_names, {})
            >> Option('l', "lefts", left_reads, {})
            >> Option('r', "rights", right_reads, {})
            >> Option('o', "out", out_root, "")
            >> Option('p', "dump-annotation", propagation_dump, "")
            >> Option('e', "dump-edges", edges_dump, "")
        ;
        if (sample_names.empty() == left_reads.empty()  && //All options of this group
            left_reads.empty()   == right_reads.empty() && //must simultaneously present or not
            right_reads.empty()  == out_root.empty())
            no_binning = sample_names.empty();
        else {
            cerr << "All of -n -l -r -o options must present for read binning" << endl;
            throw OptionNotFoundEx();
        }
    } catch(GetOptEx &ex) {
        cerr << "Usage: prop_binning -k <K> -s <saves path> -c <contigs path> -f <splits path> "
                "-a <binning annotation> [-t <length threshold>] [-b <bins to propagate>] "
                "[-n <sample names> -l <left reads> -r <right reads> -o <reads output root>] "
                "[-p <propagation info dump>] [-e <propagated edges dump>]"  << endl;
        exit(1);
    }

    vector<bin_id> bins_of_interest;
    if (!bins_file.empty()) {
        ifstream bins_stream(bins_file);
        bin_id bin;
        while (!bins_stream.eof()) {
            bins_stream >> bin;
            bins_of_interest.push_back(bin);
            bins_stream.ignore(numeric_limits<std::streamsize>::max(), '\n'); //Skip the rest of bin info
        }
        INFO("Loaded " << bins_of_interest.size() << " interesting bins");
    }

    conj_graph_pack gp(k, "tmp", 1);
    gp.kmer_mapper.Attach();

    INFO("Load graph and clustered paired info from " << saves_path);
    graphio::ScanWithClusteredIndices(saves_path, gp, gp.clustered_indices);

    //Propagation stage
    INFO("Using contigs from " << contigs_path);
    io::FileReadStream contigs_stream(contigs_path);
    io::FileReadStream split_stream(splits_path);

    AnnotationStream annotation_in(annotation_path);

    AnnotationFiller filler(gp, bins_of_interest);
    EdgeAnnotation edge_annotation = filler(contigs_stream, split_stream, annotation_in);

    INFO("Propagation launched");
    AnnotationPropagator propagator(gp, length_threshold);
    propagator.Run(contigs_stream, edge_annotation);
    INFO("Propagation finished");

    if (!edges_dump.empty()) {
        INFO("Dumping propagated edges to " << edges_dump);
        DumpEdges(gp.g, edges_dump);
    }

    if (!propagation_dump.empty()) {
        INFO("Dumping propagated annotation to " << propagation_dump);
        DumpAnnotation(gp.g, edge_annotation, propagation_dump);
    }

    //Binning stage
    if (!no_binning) {
        INFO("Binning reads into " << out_root);
        for (size_t i = 0; i < sample_names.size(); ++i)
            BinReads(gp, out_root, sample_names[i], left_reads[i], right_reads[i], edge_annotation, bins_of_interest);
    }
    return 0;
}
