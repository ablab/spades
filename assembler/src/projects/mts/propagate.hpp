//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/reads/single_read.hpp"
#include "io/reads_io/io_helper.hpp"
#include "io/reads_io/osequencestream.hpp"
#include "annotation.hpp"

namespace debruijn_graph {

inline void DumpEdgesAndAnnotation(const Graph& g,
                                   const EdgeAnnotation& edge_annotation,
                                   const string& out_prefix) {
    io::osequencestream oss(out_prefix + ".fasta");
    AnnotationOutStream annotation_out(out_prefix + ".ann");
    for (auto it = g.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        io::SingleRead edge_read("NODE_" + ToString(g.int_id(e)),
                                 g.EdgeNucls(e).str());
        oss << edge_read;
        auto relevant_bins = edge_annotation.RelevantBins(edge_read);
        if (!relevant_bins.empty()) {
            annotation_out << ContigAnnotation(GetId(edge_read),
                                               vector<bin_id>(relevant_bins.begin(), relevant_bins.end()));
        }
    }
}

class AnnotationPropagator {
    const conj_graph_pack& gp_;

public:
    AnnotationPropagator(const conj_graph_pack& gp) :
                     gp_(gp) {
    }

    void Run(io::SingleStream& contigs, EdgeAnnotation& edge_annotation);
};

}
