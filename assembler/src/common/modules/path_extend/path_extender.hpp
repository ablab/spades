//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* Copyright (c) 2014-2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "extension_chooser.hpp"
#include "overlap_analysis.hpp"
#include "path_filter.hpp"
#include "gap_analyzer.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"

#include <cmath>

namespace path_extend {

template<typename... Args>
inline BidirectionalPath& CreatePath(PathContainer &paths,
                                     GraphCoverageMap &coverage_map,
                                     Args&&... args) {
    auto p = paths.CreatePair(std::forward<Args>(args)...);
    coverage_map.Subscribe(p);

    return p.first;
}

inline BidirectionalPath& AddPath(PathContainer &paths,
                                  std::unique_ptr<BidirectionalPath> path,
                                  GraphCoverageMap &coverage_map) {
    auto p = paths.Add(std::move(path));
    coverage_map.Subscribe(p);

    return p.first;
}

class ShortLoopEstimator {
public:
    //Path must end with forward cycle edge, contain at least 2 edges and must not contain backward cycle edges
    //Returns 0 (for no loops), 1 (for a single loop) or 2 (for many loops)
    virtual size_t EstimateSimpleCycleCount(const BidirectionalPath& path, EdgeId backward_edge, EdgeId exit_edge) const = 0;

    virtual ~ShortLoopEstimator() {};
};

class ShortLoopResolver {
public:
    static const size_t BASIC_N_CNT = 100;

    ShortLoopResolver(const Graph& g, std::shared_ptr<ShortLoopEstimator> loop_estimator)
            : g_(g), loop_estimator_(loop_estimator) { }

    void ResolveShortLoop(BidirectionalPath& path) const {
        EdgeId back_cycle_edge;
        EdgeId loop_outgoing;
        EdgeId loop_incoming;
        if (path.Size() >=1 && GetLoopAndExit(g_, path.Back(), back_cycle_edge, loop_outgoing, loop_incoming)) {
            DEBUG("Resolving short loop...");
            MakeBestChoice(path, back_cycle_edge, loop_outgoing, loop_incoming);
            DEBUG("Resolving short loop done");
        }
    }

private:
    DECL_LOGGER("PathExtender")
    const Graph& g_;
    std::shared_ptr<ShortLoopEstimator> loop_estimator_;

    void UndoCycles(BidirectionalPath& p, EdgeId back_cycle_edge, EdgeId loop_incoming) const {
        if (p.Size() <= 2) {
            return;
        }
        EdgeId forward_cycle_edge = p.Back();
        size_t loop_start_index = p.Size();
        while (loop_start_index > 2) {
            if (p.At(loop_start_index - 1) == forward_cycle_edge && p.At(loop_start_index - 2) == back_cycle_edge) {
                loop_start_index -= 2;
            } else {
                break;
            }
        }

        if (p.At(loop_start_index - 1) == forward_cycle_edge) {
            p.PopBack(p.Size() - loop_start_index);
        } else if (loop_start_index != p.Size()) {
            //Means we jumped into the loop
            DEBUG("Jumped inside the loop, back loop edge " << g_.int_id(back_cycle_edge) << ", forward loop edge " << g_.int_id(forward_cycle_edge) << ", loop starts @ " << loop_start_index);
            p.PrintDEBUG();
            Gap gap = p.GapAt(loop_start_index);
            p.PopBack(p.Size() - loop_start_index);

            if (p.Back() != loop_incoming) {
                p.PushBack(loop_incoming,
                           Gap(std::max(0, gap.gap - (int)g_.length(loop_incoming) - (int) g_.length(forward_cycle_edge)),
                               {gap.trash.previous, 0}));
            }
            p.PushBack(forward_cycle_edge);
            DEBUG("Restored the path");
            p.PrintDEBUG();
        }
    }

    //edges -- first edge is loop's back edge, second is loop exit edge
    void MakeBestChoice(BidirectionalPath& path, EdgeId back_cycle_edge, EdgeId loop_outgoing, EdgeId loop_incoming) const {
        EdgeId forward_cycle_edge = path.Back();
        UndoCycles(path, back_cycle_edge, loop_incoming);

        //Expects 0 (for no loops), 1 (for a single loop) or 2 (for many loops, will insert back_cycle_edge and Ns)
        size_t loop_count = loop_estimator_->EstimateSimpleCycleCount(path, back_cycle_edge, loop_outgoing);
        if (loop_count > 0) {
            path.PushBack(back_cycle_edge);
            if (loop_count == 1) {
                DEBUG("Single loop");
                path.PushBack(forward_cycle_edge);
                path.PushBack(loop_outgoing);
            }
            else {
//TODO:: what should be here?
                if (cfg::get().pd && loop_count * (g_.length(forward_cycle_edge) + g_.length(loop_outgoing)) < 1000) {
                    INFO(" Plasmid mode: full loop resolving. Loop multiplicity: " << loop_count);
                    INFO(" Loop edges " << forward_cycle_edge << " " << back_cycle_edge);
                    for(size_t i = 0; i < loop_count - 1; i++) {
                        path.PushBack(forward_cycle_edge);
                        path.PushBack(back_cycle_edge);
                    }
                    path.PushBack(forward_cycle_edge);
                    path.PushBack(loop_outgoing);
                } else {
                    DEBUG("Multiple cycles");
                    //If the forward edge is shorter than K, avoid overlapping bases between backward edge and outgoing edge
                    //Make sure that the N-stretch will be exactly 100 bp
                    uint32_t overlapping_bases = (uint32_t) std::max(int(g_.k()) - int(g_.length(forward_cycle_edge)), 0);
                    path.PushBack(loop_outgoing, Gap(int(g_.k() + BASIC_N_CNT - overlapping_bases), {0, overlapping_bases}));

                }
            }
        }
        else {
            path.PushBack(loop_outgoing);
        }
    }
};

class CoverageLoopEstimator : public ShortLoopEstimator {
public:
    CoverageLoopEstimator(const Graph& g, const omnigraph::FlankingCoverage<Graph>& flanking_cov)
            : g_(g), flanking_cov_(flanking_cov) {

    }

    //Path must end with forward cycle edge, contain at least 2 edges and must not contain backward cycle edges
    //Returns 0 (for no loops), 1 (for a single loop) or 2 (for many loops)
    size_t EstimateSimpleCycleCount(const BidirectionalPath& path, EdgeId backward_edge, EdgeId exit_edge) const override {
        VERIFY(path.Size() > 1);
        EdgeId forward_edge = path.Back();
        EdgeId incoming_edge = path[path.Size() - 2];
        double in_cov = flanking_cov_.GetOutCov(incoming_edge);
        double out_cov = flanking_cov_.GetInCov(exit_edge);
        double avg_coverage = (in_cov + out_cov) / 2.0;

        double fwd_count = math::round(g_.coverage(forward_edge) / avg_coverage);
        double back_count = math::round(g_.coverage(backward_edge) / avg_coverage);
        size_t result = (size_t) math::round(std::max(0.0, std::min(fwd_count - 1.0, back_count)));

        DEBUG("loop with start " << g_.int_id(incoming_edge)
                <<" e1 " << g_.int_id(forward_edge)
                << " e2 " << g_.int_id(backward_edge)
                << " out " <<g_.int_id(exit_edge)
                << " cov in = " << in_cov
                << " cov out " << out_cov
                << " cov " << avg_coverage
                << " cov e1 = " << g_.coverage(forward_edge)
                << " cov e2 = " << g_.coverage(backward_edge)
                << " fwd_count = " << fwd_count
                << " back_count = " << back_count
                << " result = " << result);

        return result;
    }

private:
    const Graph& g_;
    const omnigraph::FlankingCoverage<Graph>& flanking_cov_;
};

class PairedInfoLoopEstimator: public ShortLoopEstimator {
    const Graph& g_;
    std::shared_ptr<WeightCounter> wc_;
    double weight_threshold_;

public:
    PairedInfoLoopEstimator(const Graph &g, std::shared_ptr<WeightCounter> wc, double weight_threshold = 0.0)
            : g_(g),
              wc_(wc),
              weight_threshold_(weight_threshold) { }

    //Path must end with forward cycle edge, contain at least 2 edges and must not contain backward cycle edges
    //Returns 0 (for no loops), 1 (for a single loop) or 2 (for many loops)
    size_t EstimateSimpleCycleCount(const BidirectionalPath& path, EdgeId backward_edge, EdgeId /*exit_edge*/) const override {
        VERIFY(path.Size() > 1);
        VERIFY(wc_ != nullptr);
        EdgeId forward_cycle_edge = path.Back();

        size_t result = 0;
        double lopp_edge_weight = wc_->CountWeight(path, backward_edge);
        if (math::gr(lopp_edge_weight, weight_threshold_)) {
            //Paired information on loop back edges exits => at leat one iteration
            //Looking for paired information supporting more than 1 cycle
            if (NoSelfPairedInfo(backward_edge, forward_cycle_edge)) {
                //More likely to be a single cycle
                DEBUG("Single loop");
                result = 1;
            }
            else {
                DEBUG("Multiple cycles");
                //More likely to be a 2 or more cycles
                result = 2;
            }
        }
        return result;
    }

private:
    bool NoSelfPairedInfo(EdgeId back_cycle_edge, EdgeId forward_cycle_edge) const {
        size_t is = wc_->PairedLibrary().GetISMax();
        int forward_len = (int) g_.length(forward_cycle_edge);
        bool exists_pi = true;

        // FIXME: rethink logic
        auto cycle = BidirectionalPath::create(g_, back_cycle_edge);
        while (cycle->Length() < is + g_.length(back_cycle_edge)) {
            auto w = wc_->CountWeight(*cycle, back_cycle_edge, std::set<size_t>(), forward_len);
            if (math::gr(w, weight_threshold_)) {
                //Paired information found within loop
                DEBUG("Found PI with back weight " << w << ", weight threshold " << weight_threshold_);
                exists_pi = false;
                break;
            }
            cycle->PushBack(back_cycle_edge, Gap(forward_len));
        }

        return exists_pi;
    }
};

class CombinedLoopEstimator: public ShortLoopEstimator {
public:
    CombinedLoopEstimator(const Graph& g,
                          const omnigraph::FlankingCoverage<Graph>& flanking_cov,
                          std::shared_ptr<WeightCounter> wc,
                          double weight_threshold = 0.0)
        : pi_estimator_(g, wc, weight_threshold),
          cov_estimator_(g, flanking_cov) {}

    //Path must end with forward cycle edge, contain at least 2 edges and must not contain backward cycle edges
    //Returns 0 (for no loops), 1 (for a single loop) or 2 (for many loops)
    size_t EstimateSimpleCycleCount(const BidirectionalPath& path, EdgeId backward_edge, EdgeId exit_edge) const override {
        size_t result = pi_estimator_.EstimateSimpleCycleCount(path, backward_edge, exit_edge);
        if (result == 1) {
            //Verify using coverage
            if (cov_estimator_.EstimateSimpleCycleCount(path, backward_edge, exit_edge) > 1)
                result = 2;
        }
        return result;
    }

private:
    PairedInfoLoopEstimator pi_estimator_;
    CoverageLoopEstimator cov_estimator_;
};

//Detects a cycle as a minsuffix > IS present earlier in the path. Overlap is allowed.
class InsertSizeLoopDetector {
protected:
    GraphCoverageMap visited_cycles_coverage_map_;
    PathContainer path_storage_;
    size_t min_cycle_len_;

public:
    InsertSizeLoopDetector(const Graph& g, size_t is):
        visited_cycles_coverage_map_(g),
        path_storage_(),
        min_cycle_len_(is) {
    }

    bool CheckCycledNonIS(const BidirectionalPath& path) const {
        if (path.Size() <= 2) {
            return false;
        }
        BidirectionalPath last = path.SubPath(path.Size() - 2);
        int pos = path.FindFirst(last);
        VERIFY(pos >= 0);
        return size_t(pos) != path.Size() - 2;
    }

    /// @returns whether there is more than one inclusion of the path tail that is not shorter than 'min_cycle_len'.
    bool CheckCycled(const BidirectionalPath& path) const {
        return FindCycleStart(path) != -1;
    }

    /// @returns the starting position of the first suffix that is longer or equal to 'min_cycle_len'.
    ///          Returns '-1' if the entire path is shorter than 'min_cycle_len'.
    int FindPosIS(const BidirectionalPath& path) const {
        int i = (int) path.Size() - 1;
        while (i >= 0 && path.LengthAt(i) < min_cycle_len_) {
            --i;
        }
        return i;
    }

    int FindCycleStart(const BidirectionalPath& path) const {
        TRACE("Looking for IS cycle " << min_cycle_len_);
        int i = FindPosIS(path);
        TRACE("last is pos " << i);
        if (i < 0) return -1;
//Tail
        BidirectionalPath last = path.SubPath(i);

        int pos = path.FindFirst(last);
// not cycle
        if (pos == i) pos = -1;
        TRACE("looking for 1sr IS cycle " << pos);
        return pos;
    }

//After cycle detected, removes min suffix > IS.
//returns the beginning of the cycle.
    int RemoveCycle(BidirectionalPath& path) const {
        int pos = FindCycleStart(path);
        DEBUG("Found IS cycle " << pos);
        if (pos == -1) {
            return -1;
        }

        int last_edge_pos = FindPosIS(path);
        VERIFY(last_edge_pos > -1);
        DEBUG("last edge pos " << last_edge_pos);
        VERIFY(last_edge_pos > pos);
        path.PopBack(path.Size() - (size_t)last_edge_pos);

        VERIFY((int) path.Size() == last_edge_pos);
        VERIFY(pos < (int) path.Size());
        DEBUG("result pos " <<pos);
        return pos;
    }

    //seems that it is outofdate
    bool InExistingLoop(const BidirectionalPath& path) {
        DEBUG("Checking existing loops");
        for (const auto &entry : visited_cycles_coverage_map_.GetEdgePaths(path.Back())) {
            const BidirectionalPath &cycle = *entry.first;
            DEBUG("checking  cycle ");
            int pos = path.FindLast(cycle);
            if (pos == -1)
                continue;

            int start_cycle_pos = pos + (int) cycle.Size();
            bool only_cycles_in_tail = true;
            int last_cycle_pos = start_cycle_pos;
            DEBUG("start_cycle pos "<< last_cycle_pos);
            for (int i = start_cycle_pos; i < (int) path.Size() - (int) cycle.Size(); i += (int) cycle.Size()) {
                if (!path.CompareFrom(i, cycle)) {
                    only_cycles_in_tail = false;
                    break;
                } else {
                    last_cycle_pos = i + (int) cycle.Size();
                    DEBUG("last cycle pos changed " << last_cycle_pos);
                }
            }
            DEBUG("last_cycle_pos " << last_cycle_pos);
            only_cycles_in_tail = only_cycles_in_tail && cycle.CompareFrom(0, path.SubPath(last_cycle_pos));
            if (only_cycles_in_tail) {
// seems that most of this is useless, checking
                VERIFY (last_cycle_pos == start_cycle_pos);
                DEBUG("find cycle " << last_cycle_pos);
                DEBUG("path");
                path.PrintDEBUG();
                DEBUG("last subpath");
                path.SubPath(last_cycle_pos).PrintDEBUG();
                DEBUG("cycle");
                cycle.PrintDEBUG();
                DEBUG("last_cycle_pos " << last_cycle_pos << " path size " << path.Size());
                VERIFY(last_cycle_pos <= (int)path.Size());
                DEBUG("last cycle pos + cycle " << last_cycle_pos + (int)cycle.Size());
                VERIFY(last_cycle_pos + (int)cycle.Size() >= (int)path.Size());

                return true;
            }
        }
        return false;
    }

    void AddCycledEdges(const BidirectionalPath& path, size_t pos) {
        if (pos >= path.Size()) {
            DEBUG("Wrong position in IS cycle");
            return;
        }

        auto p = path_storage_.CreatePair(path.SubPath(pos));

        visited_cycles_coverage_map_.Subscribe(p);
        DEBUG("add cycle");
        p.first.PrintDEBUG();
    }
};

class PathExtender {
public:
    explicit PathExtender(const Graph &g)
            : g_(g) { }

    virtual ~PathExtender() = default;
    virtual bool MakeGrowStep(BidirectionalPath& path, PathContainer* paths_storage = nullptr) = 0;

protected:
    const Graph &g_;
    DECL_LOGGER("PathExtender")
};


class CompositeExtender {
private:
    bool MakeGrowStep(BidirectionalPath& path, PathContainer* paths_storage);
    void GrowAllPaths(PathContainer& paths, PathContainer& result);

public:
    CompositeExtender(const Graph &g, GraphCoverageMap& cov_map,
                      UsedUniqueStorage &unique,
                      const std::vector<std::shared_ptr<PathExtender>> &pes)
            : g_(g),
              cover_map_(cov_map),
              used_storage_(unique),
              extenders_(pes) {}

    void GrowAll(PathContainer& paths, PathContainer& result);
    void GrowPath(BidirectionalPath& path, PathContainer* paths_storage) {
        while (MakeGrowStep(path, paths_storage)) { }
    }

private:
    const Graph &g_;
    GraphCoverageMap &cover_map_;
    UsedUniqueStorage &used_storage_;
    std::vector<std::shared_ptr<PathExtender>> extenders_;
};


//All Path-Extenders inherit this one
class LoopDetectingPathExtender : public PathExtender {
    const bool use_short_loop_cov_resolver_;
    ShortLoopResolver cov_loop_resolver_;

    InsertSizeLoopDetector is_detector_;
    UsedUniqueStorage &used_storage_;

protected:
    const bool investigate_short_loops_;
    const GraphCoverageMap &cov_map_;

    bool TryUseEdge(BidirectionalPath &path, EdgeId e, const Gap &gap);
    bool DetectCycle(BidirectionalPath& path);

    bool DetectCycleScaffolding(BidirectionalPath& path, EdgeId e) {
        auto temp_path = BidirectionalPath::clone(path);
        temp_path->PushBack(e);
        return is_detector_.CheckCycledNonIS(*temp_path);
    }

    virtual bool MakeSimpleGrowStep(BidirectionalPath& path, PathContainer* paths_storage = nullptr) = 0;
    virtual bool ResolveShortLoopByCov(BidirectionalPath& path);
    virtual bool ResolveShortLoopByPI(BidirectionalPath& path) = 0;
    virtual bool CanInvestigateShortLoop() const noexcept {
        return false;
    }

public:
    LoopDetectingPathExtender(const Graph &graph, const omnigraph::FlankingCoverage<Graph> &flanking_cov,
                              const GraphCoverageMap &cov_map, UsedUniqueStorage &unique,
                              bool investigate_short_loops, bool use_short_loop_cov_resolver,
                              size_t is)
            : PathExtender(graph),
              use_short_loop_cov_resolver_(use_short_loop_cov_resolver),
              cov_loop_resolver_(graph, std::make_shared<CoverageLoopEstimator>(graph, flanking_cov)),
              is_detector_(graph, is),
              used_storage_(unique),
              investigate_short_loops_(investigate_short_loops),
              cov_map_(cov_map)
            {}

    bool TryToResolveTwoLoops(BidirectionalPath& path);
    bool MakeGrowStep(BidirectionalPath& path, PathContainer* paths_storage) override;

private:
    bool ResolveShortLoop(BidirectionalPath& p) {
        if (use_short_loop_cov_resolver_) {
            return ResolveShortLoopByCov(p);
        } else {
            return ResolveShortLoopByPI(p);
        }
    }

    bool InvestigateShortLoop() const noexcept {
        return investigate_short_loops_ && (use_short_loop_cov_resolver_ || CanInvestigateShortLoop());
    }

    boost::optional<bool> TryToResolveShortLoop(BidirectionalPath& path);
protected:
    DECL_LOGGER("LoopDetectingPathExtender")
};

class SimpleExtender: public LoopDetectingPathExtender {

protected:
    std::shared_ptr<ExtensionChooser> extensionChooser_;
    ShortLoopResolver loop_resolver_;
    double weight_threshold_;

    void FindFollowingEdges(BidirectionalPath& path, ExtensionChooser::EdgeContainer * result);
public:
    SimpleExtender(const Graph &graph, const omnigraph::FlankingCoverage<Graph> &flanking_cov,
                   const GraphCoverageMap &cov_map, UsedUniqueStorage &unique,
                   std::shared_ptr<ExtensionChooser> ec,
                   bool investigate_short_loops, bool use_short_loop_cov_resolver,
                   size_t is, double weight_threshold = 0.0)
           : LoopDetectingPathExtender(graph, flanking_cov, cov_map, unique, investigate_short_loops, use_short_loop_cov_resolver, is)
           , extensionChooser_(ec)
           , loop_resolver_(graph, std::make_shared<CombinedLoopEstimator>(graph, flanking_cov, extensionChooser_->wc(), weight_threshold))
           , weight_threshold_(weight_threshold)
        {}

    SimpleExtender(const GraphPack &gp,
                   const GraphCoverageMap &cov_map,
                   UsedUniqueStorage &unique,
                   std::shared_ptr<ExtensionChooser> ec,
                   size_t is,
                   bool investigate_short_loops,
                   bool use_short_loop_cov_resolver,
                   double weight_threshold = 0.0)
            : SimpleExtender(gp.get<Graph>(), gp.get<omnigraph::FlankingCoverage<Graph>>(),
                             cov_map, unique, ec, investigate_short_loops, use_short_loop_cov_resolver,
                             is, weight_threshold)
        {}

    std::shared_ptr<ExtensionChooser> GetExtensionChooser() const { return extensionChooser_;  }
    bool CanInvestigateShortLoop() const noexcept override { return extensionChooser_->WeightCounterBased();  }
    bool ResolveShortLoopByPI(BidirectionalPath& path) override;

    bool MakeSimpleGrowStep(BidirectionalPath& path, PathContainer* paths_storage) override {
        ExtensionChooser::EdgeContainer candidates;
        return FilterCandidates(path, candidates) && AddCandidates(path, paths_storage, candidates);
    }

protected:
    virtual bool FilterCandidates(BidirectionalPath& path, ExtensionChooser::EdgeContainer& candidates);
    virtual bool AddCandidates(BidirectionalPath& path, PathContainer* /*paths_storage*/, ExtensionChooser::EdgeContainer& candidates);

private:    
    DECL_LOGGER("SimpleExtender")
};

class MultiExtender: public SimpleExtender {
public:
    using SimpleExtender::SimpleExtender;

protected:
    bool AddCandidates(BidirectionalPath& path, PathContainer* paths_storage, ExtensionChooser::EdgeContainer& candidates) override;

private:
    DECL_LOGGER("MultiExtender")

};


class ScaffoldingPathExtender: public LoopDetectingPathExtender {
    std::shared_ptr<ExtensionChooser> extension_chooser_;
    ExtensionChooser::EdgeContainer sources_;
    std::shared_ptr<GapAnalyzer> gap_analyzer_;
    bool avoid_rc_connections_;

//When check_sink_ set to false we can scaffold not only tips
    bool check_sink_;

    void InitSources();

    bool IsSink(EdgeId e) const {
        return g_.OutgoingEdgeCount(g_.EdgeEnd(e)) == 0;
    }

    Gap ConvertGapDescription(const GapDescription &gap) const {
        if (gap == GapDescription()) {
            return Gap::INVALID();
        }
        return Gap(gap.estimated_dist() + int(g_.k())
                   - int(gap.left_trim()) - int(gap.right_trim()),
                   {uint32_t(gap.left_trim()), uint32_t(gap.right_trim())}, false);
    }

protected:
    virtual bool CheckGap(const Gap &/*gap*/) const { return true; }
    bool ResolveShortLoopByCov(BidirectionalPath&) override { return false; }
    bool ResolveShortLoopByPI(BidirectionalPath&) override { return false;  }

    //TODO fix awful design with virtual CheckGap and must_overlap flag!
    bool MakeSimpleGrowStepForChooser(BidirectionalPath& path, std::shared_ptr<ExtensionChooser> ec,
                                      bool must_overlap = false);

    Gap NormalizeGap(Gap gap) const {
        VERIFY(gap != Gap::INVALID());
        if (gap.OverlapAfterTrim(g_.k()) > 0)
            gap.trash.current += gap.OverlapAfterTrim(g_.k());
        return gap;
    }

public:
    ScaffoldingPathExtender(const GraphPack &gp,
                            const GraphCoverageMap &cov_map,
                            UsedUniqueStorage &unique,
                            std::shared_ptr<ExtensionChooser> extension_chooser,
                            std::shared_ptr<GapAnalyzer> gap_analyzer,
                            size_t is,
                            bool investigate_short_loops,
                            bool avoid_rc_connections,
                            bool check_sink = true)
    :
            LoopDetectingPathExtender(gp.get<Graph>(), gp.get<omnigraph::FlankingCoverage<Graph>>(), cov_map, unique, investigate_short_loops, false, is),
            extension_chooser_(extension_chooser),
            gap_analyzer_(gap_analyzer),
            avoid_rc_connections_(avoid_rc_connections),
            check_sink_(check_sink) {
        InitSources();
    }

    bool MakeSimpleGrowStep(BidirectionalPath& path, PathContainer* /*paths_storage*/) override {
        return MakeSimpleGrowStepForChooser(path, extension_chooser_);
    }

    std::shared_ptr<ExtensionChooser> GetExtensionChooser() const {
        return extension_chooser_;
    }

private:
    DECL_LOGGER("ScaffoldingPathExtender");
};

class RNAScaffoldingPathExtender: public ScaffoldingPathExtender {
    std::shared_ptr<ExtensionChooser> strict_extension_chooser_;
    int min_overlap_;
    
protected:
    bool CheckGap(const Gap &gap) const override {
        return gap.OverlapAfterTrim(g_.k()) >= min_overlap_;
    }

public:
    RNAScaffoldingPathExtender(const GraphPack &gp,
                               const GraphCoverageMap &cov_map,
                               UsedUniqueStorage &unique,
                               std::shared_ptr<ExtensionChooser> extension_chooser,
                               std::shared_ptr<ExtensionChooser> strict_extension_chooser,
                               std::shared_ptr<GapAnalyzer> gap_joiner,
                               size_t is,
                               bool investigate_short_loops,
                               int min_overlap = 0)
    :
            ScaffoldingPathExtender(gp, cov_map, unique, extension_chooser, gap_joiner, is, investigate_short_loops, true),
            strict_extension_chooser_(strict_extension_chooser), min_overlap_(min_overlap) {}

    bool MakeSimpleGrowStep(BidirectionalPath& path, PathContainer* /*paths_storage*/) override;

private:
    DECL_LOGGER("RNAScaffoldingPathExtender");
};


}
