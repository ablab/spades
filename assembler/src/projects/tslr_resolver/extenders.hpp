#pragma once

#include "tslr_extension_chooser.hpp"

using namespace path_extend;

namespace barcode_index {
/*This class is similar to CompositeExtender, but
  with slightly different usage of coverage map.
  Huge code duplication need to be removed */

    class PathJoiner : public ContigsMaker {
    public:
        PathJoiner(const Graph & g,
                        GraphCoverageMap& cov_map,
                        size_t max_diff_len,
                        size_t max_repeat_length,
                        bool detect_repeats_online,
                        std::unordered_map <size_t, size_t> id_to_index)
            : ContigsMaker(g),
            cover_map_(cov_map),
            repeat_detector_(g, cover_map_, 2 * max_repeat_length),
            extenders_(),
            max_diff_len_(max_diff_len),
            max_repeat_len_(max_repeat_length),
            detect_repeats_online_(detect_repeats_online),
            id_to_index_(id_to_index) {}

        PathJoiner(const Graph & g,
                        GraphCoverageMap& cov_map,
                        vector<shared_ptr<PathExtender> > pes,
                        const ScaffoldingUniqueEdgeStorage& unique,
                        size_t max_diff_len,
                        size_t max_repeat_length,
                        bool detect_repeats_online,
                        std::unordered_map <size_t, size_t>& id_to_index)
                : ContigsMaker(g),
                cover_map_(cov_map),
                repeat_detector_(g, cover_map_, 2 * max_repeat_length),
                extenders_(),
                max_diff_len_(max_diff_len),
                max_repeat_len_(max_repeat_length),
                detect_repeats_online_(detect_repeats_online),
                id_to_index_(id_to_index) {
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
        unordered_map <size_t, size_t>& id_to_index_;

        void SubscribeCoverageMap(BidirectionalPath *path) {
            path->Subscribe(&cover_map_);
            for (size_t i = 0; i < path->Size(); ++i) {
                cover_map_.BackEdgeAdded(path->At(i), path, path->GapAt(i));
            }
        }

        void GrowAllPaths(PathContainer &paths, PathContainer &result) {
            size_t sum_length = 0;
            cover_map_.Clear();
            for (size_t i = 0; i < paths.size(); ++i) {
                BidirectionalPath *path = new BidirectionalPath(*paths.Get(i));
                DEBUG("new path id " << paths.Get(i) -> GetId())
                BidirectionalPath *conjugatePath = new BidirectionalPath(*paths.GetConjugate(i));
                DEBUG("conjugate path id " << paths.Get(i) -> GetId())
                result.AddPair(path, conjugatePath);
                id_to_index_[path->GetId()] = result.size() - 1;
                id_to_index_[conjugatePath->GetId()] = result.size() - 1;
                SubscribeCoverageMap(path);
                SubscribeCoverageMap(conjugatePath);
                path->CheckConjugateEnd(max_repeat_len_);
                DEBUG("result path " << path->GetId());
                sum_length += path->Length();
                DEBUG("path length " << path->Length());
                path->Print();
            }
            DEBUG("overall path length before " << sum_length);
            sum_length = 0;
            for (size_t i = 0; i < result.size(); i++) {
                if (result.Get(i)->Size() > 0) {
                    DEBUG(i);
                    DEBUG("path id " << paths.Get(i) -> GetId());
                    VERBOSE_POWER_T2(i, 100,
                                     "Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size()
                                                  << "%)");
                    if (paths.size() > 10 && i % (paths.size() / 10 + 1) == 0) {
                        INFO("Processed " << i << " paths from " << paths.size() << " (" << i * 100 / paths.size() << "%)");
                    }
                    BidirectionalPath *path = result.Get(i);
                    BidirectionalPath *conjugatePath = result.GetConjugate(i);
                    DEBUG ("Reading ids " << path->GetId() << " " << conjugatePath->GetId() << " "
                                           << id_to_index_.size());
                    size_t count_trying = 0;
                    size_t current_path_len = 0;
                    do {
                        current_path_len = path->Length();
                        count_trying++;
                        GrowPath(*path, &result);
                        GrowPath(*conjugatePath, &result);
                    } while (count_trying < 10 && (path->Length() != current_path_len));
                    path->CheckConjugateEnd(cfg::get().max_repeat_length);
                    DEBUG("result path " << path->GetId());
                    sum_length += path->Length();
                    DEBUG("path length " << path->Length() << endl);
                    path->Print();
                }
                else {
                    DEBUG("Ignoring empty path")
                }
            }
            result.FilterEmptyPaths();
        }
        DECL_LOGGER("PathJoiner")
    };

    class ReadCloudMergingExtender : public LoopDetectingPathExtender { //Traverse forward to find long edges

    protected:

        shared_ptr<ExtensionChooser> extensionChooser_;
        size_t distance_bound_;
        double barcode_threshold_;
        size_t barcode_len_;
        size_t edge_threshold_;
        size_t max_barcodes_on_edge_;
        ScaffoldingUniqueEdgeStorage unique_storage_;
        shared_ptr<barcode_index::BarcodeMapper> mapper_;
        std::unordered_map <size_t, size_t>& id_to_index_;


        void FindFollowingEdges(BidirectionalPath &path, ExtensionChooser::EdgeContainer *result) {
            result->clear();

            //find long unique edge earlier in path
            pair<EdgeId, int> last_unique =
                    ReadCloudExtensionChooser::FindLastUniqueInPath(path, unique_storage_);
            bool long_single_edge_exists = false;
            EdgeId decisive_edge;
            if (last_unique.second != -1) {
                long_single_edge_exists = true;
                decisive_edge = last_unique.first;
            }

            if (!long_single_edge_exists) {
                DEBUG("Couldn't find single long edge");
                return;
            }

            //Check if there is too much barcodes aligned to the last unique edge
//            if (mapper_->GetTailBarcodeNumber(decisive_edge) >
//                    GetMaximalBarcodeNumber(cfg::get().ts_res.read_cloud_dataset)) {
//                DEBUG("Too many barcodes mapped to the decisive edge")
//                return;
//            }

            DEBUG("At edge " << path.Back().int_id());
            DEBUG("Decisive edge " << decisive_edge.int_id());
            DEBUG("Decisive edge barcodes: " << mapper_->GetTailBarcodeNumber(decisive_edge));

            //find reliable unique edges further in graph
            vector <EdgeId> candidates;
            auto put_checker = BarcodePutChecker<Graph>(g_, edge_threshold_, barcode_threshold_, barcode_len_,
                                                        mapper_, decisive_edge, unique_storage_, candidates);
            auto dij = BarcodeDijkstra<Graph>::CreateBarcodeBoundedDijkstra(g_, distance_bound_, put_checker);
            dij.Run(g_.EdgeEnd(path.Back()));
            result->reserve(candidates.size());
            for (auto edge : candidates) {
                result->push_back(EdgeWithDistance(edge, dij.GetDistance(g_.EdgeStart(edge))));
            }
        }

        //todo should be precounted at barcode map construction stage
        size_t GetMaximalBarcodeNumber (const std::string& path_to_tslr_dataset) const {
            size_t result = 0;
            std::ifstream fin;
            fin.open(path_to_tslr_dataset);
            string line;
            while (getline(fin, line)) {
                ++result;
            }
            return result / 2;
        }


    public:

        ReadCloudMergingExtender(const conj_graph_pack &gp,
                                 const GraphCoverageMap &cov_map,
                                 shared_ptr<ExtensionChooser> ec,
                                 size_t is,
                                 size_t max_loops,
                                 bool investigate_short_loops,
                                 bool use_short_loop_cov_resolver,
                                 size_t distance_bound,
                                 double barcode_threshold,
                                 size_t barcode_len,
                                 size_t edge_threshold,
                                 size_t max_barcodes_on_edge,
                                 const ScaffoldingUniqueEdgeStorage& unique_storage,
                                 std::unordered_map <size_t, size_t>& id_to_index)
                :
                LoopDetectingPathExtender(gp, cov_map, max_loops, investigate_short_loops, use_short_loop_cov_resolver,
                                          is),
                extensionChooser_(ec),
                distance_bound_(distance_bound),
                barcode_threshold_(barcode_threshold),
                barcode_len_(barcode_len),
                edge_threshold_(edge_threshold),
                max_barcodes_on_edge_(max_barcodes_on_edge),
                unique_storage_(unique_storage),
                mapper_(gp.barcode_mapper),
                id_to_index_(id_to_index) {
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
            DEBUG("Simple grow step");
            path.Print();
            DEBUG("Path size " << path.Size())
            DEBUG("Starting at vertex " << g_.EdgeEnd(path.Back()));
            FindFollowingEdges(path, &candidates);
            DEBUG("Found " << candidates.size() << " candidates");
            DEBUG("Filtering");
            candidates = extensionChooser_->Filter(path, candidates);
            DEBUG(candidates.size() << " candidates passed");
            return true;
        }

        virtual bool AddCandidates(BidirectionalPath& path, PathContainer* paths_storage,
                                   ExtensionChooser::EdgeContainer& candidates) {
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
//            EdgeId eid = candidates.back().e_;
//            if (used_storage_->UniqueCheckEnabled()) {
//                if (used_storage_->IsUsedAndUnique(eid)) {
//                    return false;
//                } else {
//                    used_storage_->insert(eid);
//                }
//            }
            auto candidate = candidates.back();
            DEBUG("Final candidate " << candidate.e_.int_id())
            auto covering_paths = GetCoveringPaths(path, candidate.e_);
            if (covering_paths.size() > 1)
            {
                DEBUG("Too many paths, can't determine the right one")
                return false;
            }
            if (covering_paths.size() == 0)
            {
                DEBUG("Couldn't find any covering paths")
                return false;
            }
                
            DEBUG("Trying to add next path")
            return MergePathsFromPreviousStage(path, covering_paths.back(), paths_storage);

        }

    protected:
        vector <BidirectionalPath*> GetCoveringPaths(const BidirectionalPath& path, const EdgeId& e) const {
            DEBUG("Getting covering paths...")
            BidirectionalPathSet cov_paths = cov_map_.GetCoveringPaths(e);
            vector <BidirectionalPath*> result;
            for (auto it = cov_paths.begin(); it != cov_paths.end(); ++it) {
                DEBUG("Processing path...")
                BidirectionalPath* cov_path = *it;
                DEBUG("processed path id " << cov_path -> GetId())
                DEBUG("current path id " << path.GetId())
                DEBUG("conj path id " << path.GetConjPath() -> GetId())
                if (*cov_path != path && *cov_path != *(path.GetConjPath())) {
                    result.push_back(cov_path);
                }
            }
            DEBUG("Found " << result.size() << " paths")
            return result;
        }

        bool MergePathsFromPreviousStage(BidirectionalPath& current_path,
                                         BidirectionalPath* additive_path,
                                         PathContainer* paths_storage) {
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
                DEBUG("Couldn't reach the beginning of the next path");
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
                size_t current_index = id_to_index_[current_path.GetId()];
                size_t added_index = id_to_index_[additive_path -> GetId()];
                DEBUG("current index " << current_index)
                DEBUG("added index " << added_index)
                DEBUG("Clearing path " << paths_storage->Get(added_index)->GetId())
                paths_storage->Get(added_index)->Clear();

                DEBUG("Successfully added path");
                DEBUG("Path length: " << current_path.Length());
                return true;
            }
        }
        DECL_LOGGER("InconsistentExtender")

    };
} //barcode_index
