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
    auto ext = path::extension(path);
    return path.substr(0, path.length() - ext.length()) + suffix + ext;
}

void DumpEdgesAndAnnotation(const Graph& g,
                            const EdgeAnnotation& edge_annotation,
                            const string& out_edges,
                            const string& out_annotation) {
    INFO("Dumping edges to " << out_edges << "; their annotation to " << out_annotation);
    io::osequencestream oss(out_edges);
    AnnotationOutStream annotation_out(out_annotation);
    for (auto it = g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        io::SingleRead edge_read("NODE_" + ToString(g.int_id(e)),
                                 g.EdgeNucls(e).str());
        oss << edge_read;
        auto relevant_bins = edge_annotation.Annotation(e);
        if (!relevant_bins.empty()) {
            annotation_out << ContigAnnotation(GetId(edge_read),
                                               vector<bin_id>(relevant_bins.begin(), relevant_bins.end()));
        }
    }
}

int main(int argc, char** argv) {
    using namespace GetOpt;

    //TmpFolderFixture fixture("tmp");
    create_console_logger();

    size_t k;
    string saves_path, contigs_path, splits_path, annotation_path;
    vector<string> sample_names, left_reads, right_reads;
    string out_root, propagation_dump;
    vector<bin_id> bins_of_interest;
    bool no_binning;
    try {
        GetOpt_pp ops(argc, argv);
        ops.exceptions_all();
        ops >> Option('k', k)
            >> Option('s', saves_path)
            >> Option('c', contigs_path)
            >> Option('f', splits_path)
            >> Option('a', annotation_path)
            >> Option('n', sample_names)
            >> Option('l', left_reads)
            >> Option('r', right_reads)
            >> Option('o', out_root)
            >> Option('d', propagation_dump, "")
            >> Option('b', bins_of_interest, {})
            >> OptionPresent('p', no_binning);
    } catch(GetOptEx &ex) {
        cout << "Usage: prop_binning -k <K> -s <saves path> -c <contigs path> -f <splits path> "
                "-a <binning annotation> -n <sample names> -l <left reads> -r <right reads> -o <output root> "
                "[-d <propagation info dump>] [-p to disable binning] [-b <bins of interest>*]"  << endl;
        exit(1);
    }

    for (const auto& bin_id : bins_of_interest) {
        VERIFY_MSG(bin_id.find_last_of(',') == std::string::npos, "Specify bins of interest via space, not comma");
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
    AnnotationPropagator propagator(gp);
    propagator.Run(contigs_stream, edge_annotation);
    INFO("Propagation finished");

    if (!propagation_dump.empty()) {
        INFO("Dumping propagation info to " << propagation_dump);
        DumpEdgesAndAnnotation(gp.g, edge_annotation,
                               propagation_dump + ".fasta",
                               propagation_dump + ".ann");
    }

    if (no_binning) {
        INFO("Binning was disabled with -p flag");
        return 0;
    }
    //Binning stage
//    contigs_stream.reset();
//    INFO("Using propagated annotation from " << propagated_path);
//    AnnotationStream binning_stream(propagated_path);
    for (size_t i = 0; i < sample_names.size(); ++i) {
        ContigBinner binner(gp, edge_annotation, out_root, sample_names[i]);
        INFO("Initializing binner for " << sample_names[i]);
        auto paired_stream = io::PairedEasyStream(left_reads[i], right_reads[i], false, 0);
        INFO("Running binner on " << left_reads[i] << " and " << right_reads[i]);
        binner.Run(*paired_stream);
        binner.close();
    }

    return 0;
}
