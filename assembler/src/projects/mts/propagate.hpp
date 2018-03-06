//***************************************************************************
//* Copyright (c) 2015-2016 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/reads/single_read.hpp"
#include "io/reads/io_helper.hpp"
#include "io/reads/osequencestream.hpp"
#include "annotation.hpp"

namespace debruijn_graph {

class AnnotationPropagator {
    const conj_graph_pack& gp_;
    size_t length_threshold_;

public:
    AnnotationPropagator(const conj_graph_pack& gp, size_t length_threshold) :
                     gp_(gp), length_threshold_(length_threshold) {
    }

    template<typename Result, typename... Args>
    std::shared_ptr<Result> MakePropagator(Args... args) {
        return std::make_shared<Result>(gp_, length_threshold_, args...);
    }

    void Run(io::SingleStream& contigs, EdgeAnnotation& edge_annotation);

private:
    DECL_LOGGER("AnnotationPropagator");
};

}
