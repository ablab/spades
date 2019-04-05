//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "find_best_path.hpp"
#include "hmmpath.hpp"

// PathSet<DebruijnGraphCursor> find_best_path(const hmm::Fees &fees, const std::vector<DebruijnGraphCursor> &initial,
//                                             DebruijnGraphCursor::Context context) {
//     return impl::find_best_path(fees, initial, context);
// }

// PathSet<ReversedGraphCursor<DebruijnGraphCursor>> find_best_path_rev(
//     const hmm::Fees &fees, const std::vector<ReversedGraphCursor<DebruijnGraphCursor>> &initial,
//     ReversedGraphCursor<DebruijnGraphCursor>::Context context) {
//     return impl::find_best_path(fees, initial, context);
// }
//
// PathSet<AAGraphCursor<DebruijnGraphCursor>> find_best_path(
//     const hmm::Fees &fees, const std::vector<AAGraphCursor<DebruijnGraphCursor>> &initial,
//     AAGraphCursor<DebruijnGraphCursor>::Context context) {
//     return impl::find_best_path(fees, initial, context);
// }
//
// PathSet<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>> find_best_path(
//     const hmm::Fees &fees, const std::vector<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>> &initial,
//     OptimizedRestrictedGraphCursor<DebruijnGraphCursor>::Context context) {
//     return impl::find_best_path(fees, initial, context);
// }
//
// PathSet<AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>> find_best_path(
//     const hmm::Fees &fees,
//     const std::vector<AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>> &initial,
//     AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>::Context context) {
//     return impl::find_best_path(fees, initial, context);
// }

PathSet<StringCursor> find_best_path(const hmm::Fees &fees, const std::vector<StringCursor> &initial,
                                     StringCursor::Context context) {
    return impl::find_best_path(fees, initial, context);
}

PathSet<CachedCursor> find_best_path(const hmm::Fees &fees, const std::vector<CachedCursor> &initial,
                                     CachedCursor::Context context) {
    return impl::find_best_path(fees, initial, context);
}

double score_sequence(const hmm::Fees &fees, const std::string &seq) {
    StringCursor start(0), finish(seq.length() - 1);

    auto result = find_best_path(fees, {start}, &seq);
    result.pathlink_mutable()->set_finishes({finish});
    return result.best_score();
}

double score_subsequence(const hmm::Fees &fees, const std::string &seq) {
    std::vector<StringCursor> initial;
    for (size_t i = 0; i < seq.length(); ++i) {
        initial.emplace_back(i);
    }

    auto result = find_best_path(fees, initial, &seq);
    return result.best_score();
}
