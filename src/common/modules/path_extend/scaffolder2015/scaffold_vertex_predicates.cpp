//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_vertex_predicates.hpp"

namespace path_extend {
namespace scaffolder {

LengthChecker::LengthChecker(size_t length_threshold, const Graph &g)
    : length_threshold_(length_threshold), g_(g) {}
bool LengthChecker::Check(const ScaffoldVertex &vertex) const {
    return vertex.GetLengthFromGraph(g_) < length_threshold_;
}

AndPredicate::AndPredicate(std::shared_ptr<ScaffoldVertexPredicate> first,
                           std::shared_ptr<ScaffoldVertexPredicate> second) :
    first_(first), second_(second) {}
bool AndPredicate::Check(const ScaffoldVertex &scaffold_vertex) const {
    return first_->Check(scaffold_vertex) and second_->Check(scaffold_vertex);
}

}
}