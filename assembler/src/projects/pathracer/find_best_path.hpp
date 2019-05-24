//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pathtree.hpp"
#include "restricted_cursor.hpp"
#include "aa_cursor.hpp"
#include "string_cursor.hpp"
#include "cached_cursor.hpp"
#include "debruijn_graph_cursor.hpp"
#include "reversed_cursor.hpp"

double score_sequence(const hmm::Fees &fees, const std::string &seq);
double score_subsequence(const hmm::Fees &fees, const std::string &seq);

namespace hmm {
struct Fees;
};

PathSet<DebruijnGraphCursor> find_best_path(const hmm::Fees &fees, const std::vector<DebruijnGraphCursor> &initial,
                                            DebruijnGraphCursor::Context context);

PathSet<ReversedGraphCursor<DebruijnGraphCursor>> find_best_path_rev(
    const hmm::Fees &fees, const std::vector<ReversedGraphCursor<DebruijnGraphCursor>> &initial,
    ReversedGraphCursor<DebruijnGraphCursor>::Context context);

PathSet<AAGraphCursor<DebruijnGraphCursor>> find_best_path(
    const hmm::Fees &fees, const std::vector<AAGraphCursor<DebruijnGraphCursor>> &initial,
    AAGraphCursor<DebruijnGraphCursor>::Context context);

PathSet<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>> find_best_path(
    const hmm::Fees &fees, const std::vector<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>> &initial,
    OptimizedRestrictedGraphCursor<DebruijnGraphCursor>::Context context);

PathSet<AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>> find_best_path(
    const hmm::Fees &fees,
    const std::vector<AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>> &initial,
    AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>::Context context);

PathSet<StringCursor> find_best_path(const hmm::Fees &fees, const std::vector<StringCursor> &initial,
                                     StringCursor::Context context);

PathSet<CachedCursor> find_best_path(const hmm::Fees &fees, const std::vector<CachedCursor> &initial,
                                     CachedCursor::Context context);
PathSet<AAGraphCursor<StringCursor>> find_best_path(const hmm::Fees &fees, const std::vector<AAGraphCursor<StringCursor>> &initial,
                                                    AAGraphCursor<StringCursor>::Context context);
