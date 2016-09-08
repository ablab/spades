#pragma once

#include "tslr_extension_chooser.hpp"

using namespace path_extend;

namespace tslr_resolver {
/*todo this class is similar to CompositeExtender, but
with slightly different usage of coverage map. Huge code duplication need to be removed /*/

    class PathJoiner : public ContigsMaker {
    public:
        PathJoiner(Graph & g, GraphCoverageMap& cov_map,
                size_t max_diff_len,
        size_t max_repeat_length,
        bool detect_repeats_online)
        : ContigsMaker(g),
                cover_map_(cov_map),
                repeat_detector_(g, cover_map_, 2 * max_repeat_length),
        extenders_(),
        max_diff_len_(max_diff_len),
        max_repeat_len_(max_repeat_length),
        detect_repeats_online_(detect_repeats_online) {
        }

        PathJoiner(Graph & g, GraphCoverageMap& cov_map,
                vector<shared_ptr<PathExtender> > pes,
        const ScaffoldingUniqueEdgeStorage& unique,
                size_t max_diff_len,
        size_t max_repeat_length,
        bool detect_repeats_online)
        : ContigsMaker(g),
                cover_map_(cov_map),
                repeat_detector_(g, cover_map_, 2 * max_repeat_length),
        extenders_(),
        max_diff_len_(max_diff_len),
        max_repeat_len_(max_repeat_length),
        detect_repeats_online_(detect_repeats_online) {
            extenders_ = pes;
            used_storage_ = make_shared<UsedUniqueStorage>(UsedUniqueStorage(unique));
            for (auto ex: extenders_) {
                ex->AddUniqueEdgeStorage(used_storage_);
            }
        }

        void AddExtender(shared_ptr <PathExtender> pe) {
            extenders_.push_back(pe);
            pe->AddUniqueEdgeStorage(used_storage_);
        }

        void GrowAll(PathContainer &paths, PathContainer &result) override {
            result.clear();
            GrowAllPaths(paths, result);
            LengthPathFilter filter(g_, 0);
            filter.filter(result);
        }

        void GrowPath(BidirectionalPath &path, PathContainer *paths_storage) override {
            while (MakeGrowStep(path, paths_storage)) {}
        }

        void GrowPathSimple(BidirectionalPath &path, PathContainer *paths_storage) override {
            while (MakeGrowStep(path, paths_storage)) {}
        }

        bool MakeGrowStep(BidirectionalPath &path, PathContainer *paths_storage) {
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

    private:
        GraphCoverageMap &cover_map_;
        RepeatDetector repeat_detector_;
        vector <shared_ptr<PathExtender>> extenders_;
        size_t max_diff_len_;
        size_t max_repeat_len_;
        bool detect_repeats_online_;
        shared_ptr <UsedUniqueStorage> used_storage_;

        void SubscribeCoverageMap(BidirectionalPath *path) {
            path->Subscribe(&cover_map_);
            for (size_t i = 0; i < path->Size(); ++i) {
                cover_map_.BackEdgeAdded(path->At(i), path, path->GapAt(i));
            }
        }

        void GrowAllPaths(PathContainer &paths, PathContainer &result) {
            for (size_t i = 0; i < paths.size(); ++i) {
                DEBUG(i)
                DEBUG("path id " << paths.Get(i) -> GetId())
                VERBOSE_POWER_T2(i, 100,
                                 "Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size()
                                              << "%)");
                if (paths.size() > 10 && i % (paths.size() / 10 + 1) == 0) {
                    INFO("Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size() << "%)");
                }
                if (paths.Get(i) -> Size() != 0) {

//In 2015 modes do not use a seed already used in paths.
                    if (used_storage_->UniqueCheckEnabled()) {
                        bool was_used = false;
                        for (size_t ind = 0; ind < paths.Get(i)->Size(); ind++) {
                            EdgeId eid = paths.Get(i)->At(ind);
                            if (used_storage_->IsUsedAndUnique(eid)) {
                                was_used = true;
                                break;
                            } else {
                                used_storage_->insert(eid);
                            }
                        }
                        if (was_used) {
                            DEBUG("skipping already used seed");
                            continue;
                        }
                    }
                    BidirectionalPath *path = new BidirectionalPath(*paths.Get(i));
                    DEBUG("new path id " << paths.Get(i) -> GetId())
                    BidirectionalPath *conjugatePath = new BidirectionalPath(*paths.GetConjugate(i));
                    DEBUG("conjugate path id " << paths.Get(i) -> GetId())
                    result.AddPair(path, conjugatePath);
                    SubscribeCoverageMap(path);
                    SubscribeCoverageMap(conjugatePath);
                    size_t count_trying = 0;
                    size_t current_path_len = 0;
                    do {
                        current_path_len = path->Length();
                        count_trying++;
                        GrowPath(*path, &result);
                        GrowPath(*conjugatePath, &result);
                    } while (count_trying < 10 && (path->Length() != current_path_len));
                    path->CheckConjugateEnd(max_repeat_len_);
                    DEBUG("result path " << path->GetId());
                    path->Print();
                }
            }
        }
        DECL_LOGGER("PathJoiner")
    };

    class InconsistentTSLRExtender : public LoopDetectingPathExtender { //Traverse forward to find long edges

    protected:
        shared_ptr<ExtensionChooser> extensionChooser_;
        size_t distance_bound_;
        double barcode_threshold_;
        size_t edge_threshold_;
        bool join_paths;
        ScaffoldingUniqueEdgeStorage unique_storage_;
        shared_ptr<tslr_resolver::BarcodeMapper> mapper_;
        double average_coverage_;


        void FindFollowingEdges(BidirectionalPath &path, ExtensionChooser::EdgeContainer *result) {
            result->clear();
            if (g_.OutgoingEdgeCount(g_.EdgeEnd(path.Back())) == 1) {
                result->push_back(EdgeWithDistance(*(g_.OutgoingEdges(g_.EdgeEnd(path.Back())).begin()), 0));
                return;
            }

            //find long unique edge earlier in path
            bool long_single_edge_exists = false;
            EdgeId decisive_edge;
            for (int i = static_cast<int> (path.Size()) - 1; !long_single_edge_exists && i >= 0; --i) {
                EdgeId current_edge = path.At(i);
                if (unique_storage_.IsUnique(current_edge) &&
                    g_.length(current_edge) > edge_threshold_) {
                    long_single_edge_exists = true;
                    decisive_edge = current_edge;
                }
            }

            if (!long_single_edge_exists) {
                DEBUG("Couldn't find single long edge");
                return;
            }

            if (mapper_-> GetSizeTails(decisive_edge) > 80) {
                DEBUG("Too many barcodes mapped to the decisive edge")
                return;
            }

            if (mapper_-> GetSizeTails(decisive_edge) < 5) {
                DEBUG("Not enough barcodes mapped to the decisive edge")
                return;
            }

            DEBUG("At edge " << path.Back().int_id());
            DEBUG("decisive edge " << decisive_edge.int_id());

            DEBUG("decisive edge barcodes: " <<  mapper_->GetSizeTails(decisive_edge));

            //find reliable unique edges further in graph
            vector <EdgeId> candidates;
            auto put_checker = BarcodePutChecker<Graph>(g_, edge_threshold_, barcode_threshold_,
                                                        mapper_, decisive_edge, unique_storage_, candidates);
            auto dij = BarcodeDijkstra<Graph>::CreateBarcodeBoundedDijkstra(g_, distance_bound_, put_checker);
            dij.Run(g_.EdgeEnd(path.Back()));
            result->reserve(candidates.size());
            for (auto edge : candidates) {
                result->push_back(EdgeWithDistance(edge, dij.GetDistance(g_.EdgeStart(edge))));
            }
        }


    public:

        InconsistentTSLRExtender(const conj_graph_pack &gp,
                                 const GraphCoverageMap &cov_map,
                                 shared_ptr<ExtensionChooser> ec,
                                 size_t is,
                                 size_t max_loops,
                                 bool investigate_short_loops,
                                 bool use_short_loop_cov_resolver,
                                 size_t distance_bound,
                                 double barcode_threshold,
                                 size_t edge_threshold,
                                 bool join_paths,
                                 const ScaffoldingUniqueEdgeStorage& unique_storage,
                                 double average_coverage)
                :
                LoopDetectingPathExtender(gp, cov_map, max_loops, investigate_short_loops, use_short_loop_cov_resolver,
                                          is),
                extensionChooser_(ec),
                distance_bound_(distance_bound),
                barcode_threshold_(barcode_threshold),
                edge_threshold_(edge_threshold),
                join_paths(join_paths),
                unique_storage_(unique_storage),
                mapper_(gp.barcode_mapper),
                average_coverage_(average_coverage) {
        }

        std::shared_ptr<ExtensionChooser> GetExtensionChooser() const {
            return extensionChooser_;
        }

        bool CanInvestigateShortLoop() const override {
            return extensionChooser_->WeightCounterBased();
        }

        bool ResolveShortLoopByCov(BidirectionalPath &path) override {
            return path.Size() < 1;
        }

        bool ResolveShortLoopByPI(BidirectionalPath &path) override {
            return path.Size() < 1;
        }

        bool MakeSimpleGrowStep(BidirectionalPath &path, PathContainer *paths_storage) override {
            ExtensionChooser::EdgeContainer candidates;
            return FilterCandidates(path, candidates) and AddCandidates(path, paths_storage, candidates);
        }

    protected:
        virtual bool FilterCandidates(BidirectionalPath &path, ExtensionChooser::EdgeContainer &candidates) {
            DEBUG("Cov map size " << cov_map_.size())
            if (path.Size() == 0) {
                return false;
            }
            DEBUG("Simple grow step");
            path.Print();
            DEBUG("Path size " << path.Size())
            DEBUG("Starting at vertex " << g_.EdgeEnd(path.Back()));
            FindFollowingEdges(path, &candidates);
            DEBUG("found candidates");
            DEBUG(candidates.size())
            if (candidates.size() == 1) {
                LoopDetector loop_detector(&path, cov_map_);
                if (!investigate_short_loops_ &&
                    (loop_detector.EdgeInShortLoop(path.Back()) or loop_detector.EdgeInShortLoop(candidates.back().e_))
                    && extensionChooser_->WeightCounterBased()) {
                    return false;
                }
            }
            DEBUG("more filtering");
            candidates = extensionChooser_->Filter(path, candidates);
            DEBUG("filtered candidates");
            DEBUG(candidates.size())
            return true;
        }

        virtual bool AddCandidates(BidirectionalPath &path, PathContainer *paths_storage,
                                   ExtensionChooser::EdgeContainer &candidates) {
            if (candidates.size() != 1) {
                if (candidates.size() > 1) {
                    DEBUG("Too many candidates, false");
                }
                if (candidates.size() == 0) {
                    DEBUG("No candidates found, false");
                }
                DEBUG("Final(?) path length: " << path.Length());
                return false;
            }

            DEBUG("push");
            EdgeId eid = candidates.back().e_;

            if (used_storage_->UniqueCheckEnabled()) {
                if (used_storage_->IsUsedAndUnique(eid)) {
                    return false;
                } else {
                    used_storage_->insert(eid);
                }
            }
            auto candidate = candidates.back();
            auto covering_paths = GetCoveringPaths(path, candidate.e_);
            if (!join_paths || covering_paths.size() != 1)
            {
                if (join_paths) {
                    DEBUG("Next path indetermined")
                }
                path.PushBack(eid, candidates.back().d_);
            }
                
            else {
                DEBUG("Trying to add next path")
                return AddPathFromPreviousStage(path, covering_paths.back());
            }
            
            DEBUG("Push done, true");
            DEBUG("Path length: " << path.Length());
            return true;
        }

    protected:
        vector <BidirectionalPath*> GetCoveringPaths(const BidirectionalPath& path, const EdgeId& e) const {
            DEBUG("Getting covering paths...")
            DEBUG(cov_map_.size())
            BidirectionalPathSet cov_paths = cov_map_.GetCoveringPaths(e);
            vector <BidirectionalPath*> result;
            DEBUG(cov_paths.size())
            for (auto it = cov_paths.begin(); it != cov_paths.end(); ++it) {
                DEBUG("Processing path...")
                BidirectionalPath* cov_path = *it;
                DEBUG("processed path id " << cov_path -> GetId())
                if (*cov_path != path && *cov_path != *(path.GetConjPath())) {
                    result.push_back(cov_path);
                }
            }
            return result;
        }

        bool AddPathFromPreviousStage(BidirectionalPath &current_path,
                                      BidirectionalPath* additive_path) {
            //Count distance between paths
            size_t path_len_bound = cfg::get().ts_res.topsort_bound;
            size_t overlap_size = current_path.OverlapEndSize(additive_path);
            DEBUG("First path size " << current_path.Size())
            DEBUG("Second path size " << additive_path -> Size())
            DEBUG("Overlap size " << overlap_size)
            BidirectionalPath suffix = additive_path->SubPath(overlap_size);

            VertexId start_vertex = g_.EdgeEnd(current_path.Back());
            auto dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, path_len_bound);
            dijkstra.Run(start_vertex);

            if (start_vertex != g_.EdgeStart(suffix.Front()) &&
                    !dijkstra.DistanceCounted(g_.EdgeStart(suffix.Front()))) {
                DEBUG("Couldn't reach next path")
                return false;
            }
            else {
                DEBUG("Adding path")
                DEBUG("path id " << current_path.GetId())
                DEBUG("Added path id " << additive_path -> GetId())
                DEBUG("Distance " << dijkstra.GetDistance(g_.EdgeStart(suffix.Front())))
                int gap = static_cast<int> (dijkstra.GetDistance(g_.EdgeStart(suffix.Front())));

                current_path.PushBack(suffix.At(0), gap);
                current_path.PushBack(suffix.SubPath(1));
                

                return true;
            }
        }
        DECL_LOGGER("InconsistentExtender")

    };
} //tslr_resolver
