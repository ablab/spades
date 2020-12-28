//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* Copyright (c) 2014-2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "path_extender.hpp"

namespace path_extend {

void CompositeExtender::GrowAll(PathContainer& paths, PathContainer& result) {
    result.clear();
    GrowAllPaths(paths, result);
    result.FilterEmptyPaths();
}

bool CompositeExtender::MakeGrowStep(BidirectionalPath& path, PathContainer* paths_storage) {
    DEBUG("make grow step composite extender");

    size_t current = 0;
    while (current < extenders_.size()) {
        DEBUG("step " << current << " of total " << extenders_.size());
        if (extenders_[current]->MakeGrowStep(path, paths_storage)) {
            return true;
        }
        ++current;
    }
    return false;
}

void CompositeExtender::GrowAllPaths(PathContainer& paths, PathContainer& result) {
    for (size_t i = 0; i < paths.size(); ++i) {
        VERBOSE_POWER_T2(i, 100, "Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size() << "%)");
        if (paths.size() > 10 && i % (paths.size() / 10 + 1) == 0) {
            INFO("Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size() << "%)");
        }
        //In 2015 modes do not use a seed already used in paths.
        //FIXME what is the logic here?
        if (used_storage_.UniqueCheckEnabled()) {
            bool was_used = false;
            const BidirectionalPath &p = paths.Get(i);
            for (size_t ind =0; ind < p.Size(); ind++) {
                EdgeId eid = p.At(ind);
                auto path_id = p.GetId();
                if (used_storage_.IsUsedAndUnique(eid, path_id)) {
                    DEBUG("Used edge " << g_.int_id(eid));
                    was_used = true;
                    break;
                } else {
                    used_storage_.insert(eid, path_id);
                }
            }
            if (was_used) {
                DEBUG("skipping already used seed");
                continue;
            }
        }

        if (!cover_map_.IsCovered(paths.Get(i))) {
            BidirectionalPath &path = CreatePath(result, cover_map_,
                                                 paths.Get(i));

            size_t count_trying = 0;
            size_t current_path_len = 0;
            do {
                current_path_len = path.Length();
                count_trying++;
                GrowPath(path, &result);
                GrowPath(*path.GetConjPath(), &result);
            } while (count_trying < 10 && (path.Length() != current_path_len));
                DEBUG("result path " << path.GetId());
                path.PrintDEBUG();
        }
    }
}

bool LoopDetectingPathExtender::TryUseEdge(BidirectionalPath &path, EdgeId e, const Gap &gap) {
    bool success = used_storage_.TryUseEdge(path, e, gap);
    if (success) {
        DEBUG("Adding edge. PathId: " << path.GetId() << " path length: " << path.Length() - 1 << ", fixed gap : "
              << gap.gap << ", trash length: " << gap.trash.previous << "-" << gap.trash.current);
    }
    return success;
}

bool LoopDetectingPathExtender::DetectCycle(BidirectionalPath& path) {
    DEBUG("detect cycle");
    if (is_detector_.CheckCycled(path)) {
        DEBUG("Checking IS cycle");
        int loop_pos = is_detector_.RemoveCycle(path);
        DEBUG("Removed IS cycle");
        if (loop_pos != -1) {
            is_detector_.AddCycledEdges(path, loop_pos);
            return true;
        }
    }
    return false;
}

bool LoopDetectingPathExtender::ResolveShortLoopByCov(BidirectionalPath& path) {
    if (TryToResolveTwoLoops(path))
        return true;

    LoopDetector loop_detector(path, cov_map_);
    size_t init_len = path.Length();
    bool result = false;
    while (path.Size() >= 1 && loop_detector.EdgeInShortLoop(path.Back())) {
        cov_loop_resolver_.ResolveShortLoop(path);
        if (init_len == path.Length()) {
            return result;
        } else {
            result = true;
        }
        init_len = path.Length();
    }
    return true;
}

bool LoopDetectingPathExtender::TryToResolveTwoLoops(BidirectionalPath& path) {
    EdgeId last_edge = path.Back();
    VertexId last_vertex = g_.EdgeEnd(last_edge);
    VertexId first_vertex = g_.EdgeStart(last_edge);
    DEBUG("Looking for two loop structure");
    auto between = g_.GetEdgesBetween(last_vertex, first_vertex);
    if (between.size() != 2 || g_.OutgoingEdgeCount(last_vertex) != 2 || g_.IncomingEdgeCount(first_vertex) != 2
        || g_.IncomingEdgeCount(last_vertex) != 1 || g_.OutgoingEdgeCount(first_vertex) != 1 ) {
        return false;
    }
    DEBUG("two glued cycles!");
    if (path.Size() >= 3 && path[path.Size()-3] == last_edge) {
        if (path.Size() >= 5 && path.Front() == last_edge){
            path.GetConjPath()->PopBack();
            DEBUG("removing extra");
            path.PrintDEBUG();
        }
        DEBUG("already traversed");
        return false;
    }
    //FIXME: constants
    if (g_.coverage(between[0]) > 1.3 *g_.coverage(between[1]) || g_.coverage(between[1]) > 1.3 *g_.coverage(between[0])) {
        DEBUG("coverage not close");
        return false;
    }
    if (path.Size() == 1) {
        path.PushBack(between[0]);
        path.PushBack(last_edge);
        path.PushBack(between[1]);
    } else {
        DEBUG("Resolved two edge-cycle, adding two edges");
        path.PushBack(between[0] == path[path.Size() - 2]? between[1] : between[0]);
        path.PushBack(last_edge);
    }
    return true;
}

bool LoopDetectingPathExtender::MakeGrowStep(BidirectionalPath& path, PathContainer* paths_storage) {
    if (path.IsCycle() || is_detector_.InExistingLoop(path) || DetectCycle(path))
        return false;

    if (TryToResolveTwoLoops(path))
        return true;

    if (use_short_loop_cov_resolver_) {
        auto attempt = TryToResolveShortLoop(path);
        if (attempt)
            return *attempt;
    }

    DEBUG("Making step");
    bool path_is_growed = MakeSimpleGrowStep(path, paths_storage);
    DEBUG("Made step");

    if (DetectCycle(path))
        return false;

    if (auto attempt = TryToResolveShortLoop(path))
        return *attempt;

    return path_is_growed;
}

boost::optional<bool> LoopDetectingPathExtender::TryToResolveShortLoop(BidirectionalPath& path) {
    LoopDetector loop_detector(path, cov_map_);
    if (!path.Empty() && InvestigateShortLoop()) {
        if (loop_detector.EdgeInShortLoop(path.Back())) {
            DEBUG("Edge in short loop");
            return ResolveShortLoop(path);
        }
        if (loop_detector.PrevEdgeInShortLoop()) {
            DEBUG("Prev edge in short loop");
            path.PopBack();
            return ResolveShortLoop(path);
        }
    }
    return {};
};

void SimpleExtender::FindFollowingEdges(BidirectionalPath& path, ExtensionChooser::EdgeContainer * result) {
    DEBUG("Looking for the following edges");
    result->clear();
    std::vector<EdgeId> edges;
    DEBUG("Pushing back");
    utils::push_back_all(edges, g_.OutgoingEdges(g_.EdgeEnd(path.Back())));
    result->reserve(edges.size());
    for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
        DEBUG("Adding edge w distance " << g_.int_id(*iter));
        result->push_back(EdgeWithDistance(*iter, 0));
    }
    DEBUG("Following edges found");
}


bool SimpleExtender::ResolveShortLoopByPI(BidirectionalPath& path) {
    if (extensionChooser_->WeightCounterBased()) {
        LoopDetector loop_detector(path, cov_map_);
        size_t init_len = path.Length();
        bool result = false;
        while (path.Size() >= 1 && loop_detector.EdgeInShortLoop(path.Back())) {
            loop_resolver_.ResolveShortLoop(path);
            if (init_len == path.Length()) {
                return result;
            } else {
                result = true;
            }
            init_len = path.Length();
        }
        return true;
    }
    return false;
}

bool SimpleExtender::FilterCandidates(BidirectionalPath& path, ExtensionChooser::EdgeContainer& candidates) {
    if (path.Size() == 0)
        return false;

    DEBUG("Simple grow step");
    path.PrintDEBUG();
    FindFollowingEdges(path, &candidates);
    DEBUG("found candidates");
    DEBUG(candidates.size());

    if (candidates.size() == 1) {
        LoopDetector loop_detector(path, cov_map_);
        if (!investigate_short_loops_ &&
            (loop_detector.EdgeInShortLoop(path.Back()) || loop_detector.EdgeInShortLoop(candidates.back().e_)) &&
            extensionChooser_->WeightCounterBased()) {
            return false;
        }
    }
    DEBUG("more filtering");
    candidates = extensionChooser_->Filter(path, candidates);
    DEBUG("filtered candidates");
    DEBUG(candidates.size());
    return true;
}

bool SimpleExtender::AddCandidates(BidirectionalPath& path, PathContainer* /*paths_storage*/, ExtensionChooser::EdgeContainer& candidates) {
    if (candidates.size() != 1)
        return false;

    LoopDetector loop_detector(path, cov_map_);
    DEBUG("loop detecor");
    if (!investigate_short_loops_ &&
        (loop_detector.EdgeInShortLoop(path.Back()) or loop_detector.EdgeInShortLoop(candidates.back().e_)) &&
        extensionChooser_->WeightCounterBased()) {
        return false;
    }
    DEBUG("push");
    EdgeId eid = candidates.back().e_;
    //In 2015 modes when trying to use already used unique edge, it is not added and path growing stops.
    //That allows us to avoid overlap removal hacks used earlier.
    Gap gap(std::move(candidates.back().gap_sequence_), candidates.back().d_);
    return TryUseEdge(path, eid, gap);
}

bool MultiExtender::AddCandidates(BidirectionalPath& path, PathContainer* paths_storage, ExtensionChooser::EdgeContainer& candidates) {
    if (candidates.size() == 0)
        return false;

    bool res = false;
    LoopDetector loop_detector(path, cov_map_);
    DEBUG("loop detecor");
    if (!investigate_short_loops_ &&
        (loop_detector.EdgeInShortLoop(path.Back()) or loop_detector.EdgeInShortLoop(candidates.back().e_))
        && extensionChooser_->WeightCounterBased()) {
        DEBUG("loop deteced");
        return false;
    }
    if (candidates.size() == 1) {
        DEBUG("push");
        EdgeId eid = candidates.back().e_;
        path.PushBack(eid, Gap(candidates.back().d_));
        DEBUG("push done");
        return true;
    } else if (candidates.size() == 2) {
        //Check for bulge
        auto v = g_.EdgeStart(candidates.front().e_);
        auto u = g_.EdgeEnd(candidates.front().e_);
        for (const auto& edge : candidates) {
            if (v != g_.EdgeStart(edge.e_) || u != g_.EdgeEnd(edge.e_))
                return false;
        }

        // Creating new paths for other than new candidate.
        for (size_t i = 1; i < candidates.size(); ++i) {
            DEBUG("push other candidates " << i);
            auto &p = paths_storage->Create(path);
            p.PushBack(candidates[i].e_, Gap(candidates[i].d_));
        }

        DEBUG("push");
        path.PushBack(candidates.front().e_, Gap(candidates.front().d_));
        DEBUG("push done");
        res = true;

        if (candidates.size() > 1) {
            DEBUG("Found " << candidates.size() << " candidates");
        }
    }

    return res;
}


void ScaffoldingPathExtender::InitSources() {
    sources_.clear();

    for (EdgeId e : g_.edges()) {
        if (g_.IncomingEdgeCount(g_.EdgeStart(e)) > 0)
            continue;

        sources_.push_back(EdgeWithDistance(e, 0));
    }
}

bool ScaffoldingPathExtender::MakeSimpleGrowStepForChooser(BidirectionalPath& path,
                                                           std::shared_ptr<ExtensionChooser> ec,
                                                           bool must_overlap) {
    if (path.Size() < 1 || (check_sink_ && !IsSink(path.Back()))) {
        return false;
    }

    DEBUG("Simple grow step, growing path");
    path.PrintDEBUG();
    ExtensionChooser::EdgeContainer candidates = ec->Filter(path, sources_);
    DEBUG("scaffolding candidates " << candidates.size() << " from sources " << sources_.size());

    DEBUG("Candidate size = " << candidates.size());
    if (candidates.size() != 1) {
        DEBUG("scaffolding end");
        return false;
    }

    EdgeId e = candidates.back().e_;
    if (e == path.Back() ||
        (avoid_rc_connections_ && e == g_.conjugate(path.Back()))) {
        return false;
    }

    if (this->DetectCycleScaffolding(path, e)) {
        return false;
    }

    Gap gap;
    //TODO is it ok that we either force joining or ignore its possibility
    if (check_sink_) {
        gap = ConvertGapDescription(gap_analyzer_->FixGap(GapDescription(path.Back(), e,
                                                                         candidates.back().d_ -
                                                                         int(g_.k()))));

        if (gap == Gap::INVALID()) {
            DEBUG("Looks like wrong scaffolding. PathId: "
                  << path.GetId() << " path length: " << path.Length()
                  << ", estimated gap length: " << candidates.back().d_);
            return false;
        }

        DEBUG("Gap after fixing " << gap.gap << " (was " << candidates.back().d_ << ")");
        if (must_overlap && !CheckGap(gap)) {
            DEBUG("Overlap is not large enough");
            return false;
        }
    } else {
        DEBUG("Gap joiners off");
        VERIFY(candidates.back().d_ > int(g_.k()));
        gap = Gap(candidates.back().d_, false);
    }

    return TryUseEdge(path, e, NormalizeGap(gap));
}


bool RNAScaffoldingPathExtender::MakeSimpleGrowStep(BidirectionalPath& path, PathContainer* /*paths_storage*/) {
    DEBUG("==== RNA scaffolding ===");
    path.PrintDEBUG();
    bool res = MakeSimpleGrowStepForChooser(path, GetExtensionChooser(), true);
    if (!res) {
        DEBUG("==== Second strategy ====");
        res = MakeSimpleGrowStepForChooser(path, strict_extension_chooser_);
    }
    DEBUG("==== DONE ====");
    return res;
}

}
