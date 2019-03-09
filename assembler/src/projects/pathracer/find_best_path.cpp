//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "find_best_path.hpp"
#include "cursor.hpp"
#include "hmmpath.hpp"

// PathSet<DebruijnGraphCursor> find_best_path(const hmm::Fees &fees, const std::vector<DebruijnGraphCursor> &initial,
//                                             DebruijnGraphCursor::Context context) {
//     return impl::find_best_path(fees, initial, context);
// }

// PathSet<RestrictedGraphCursor<DebruijnGraphCursor>> find_best_path(
//     const hmm::Fees &fees, const std::vector<RestrictedGraphCursor<DebruijnGraphCursor>> &initial,
//     RestrictedGraphCursor<DebruijnGraphCursor>::Context context) {
//     return impl::find_best_path(fees, initial, context);
// }
//
// PathSet<ReversalGraphCursor<DebruijnGraphCursor>> find_best_path_rev(
//     const hmm::Fees &fees, const std::vector<ReversalGraphCursor<DebruijnGraphCursor>> &initial,
//     ReversalGraphCursor<DebruijnGraphCursor>::Context context) {
//     return impl::find_best_path(fees, initial, context);
// }
//
// PathSet<AAGraphCursor<DebruijnGraphCursor>> find_best_path(
//     const hmm::Fees &fees, const std::vector<AAGraphCursor<DebruijnGraphCursor>> &initial,
//     AAGraphCursor<DebruijnGraphCursor>::Context context) {
//     return impl::find_best_path(fees, initial, context);
// }
//
// PathSet<AAGraphCursor<RestrictedGraphCursor<DebruijnGraphCursor>>> find_best_path(
//     const hmm::Fees &fees, const std::vector<AAGraphCursor<RestrictedGraphCursor<DebruijnGraphCursor>>> &initial,
//     AAGraphCursor<RestrictedGraphCursor<DebruijnGraphCursor>>::Context context) {
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
