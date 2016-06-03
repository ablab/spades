//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/reads_io/io_helper.hpp"
#include "annotation.hpp"

namespace debruijn_graph {

class AnnotationPropagator {
    const conj_graph_pack& gp_;

    void DumpContigAnnotation(io::SingleStream& contigs,
                              const EdgeAnnotation& annotation,
                              const string& annotation_out_fn) const;

public:
    AnnotationPropagator(const conj_graph_pack& gp) :
                     gp_(gp) {
    }

    void Run(io::SingleStream& contigs, EdgeAnnotation& edge_annotation);
};

}
