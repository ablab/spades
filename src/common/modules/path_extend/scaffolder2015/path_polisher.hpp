//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/basic_graph_stats.hpp"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "assembly_graph/paths/path_processor.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "modules/path_extend/paired_library.hpp"
#include "modules/path_extend/path_extender.hpp"

namespace path_extend {

class PathGapCloser {
protected:
    const Graph& g_;
    const size_t max_path_len_;
    const int min_gap_;

    virtual Gap CloseGap(const BidirectionalPath &original_path, size_t position,
                         BidirectionalPath &path) const = 0;
    DECL_LOGGER("PathGapCloser")
public:
    std::unique_ptr<BidirectionalPath> CloseGaps(const BidirectionalPath &path) const;

    PathGapCloser(const Graph& g, size_t max_path_len):
                  g_(g),
                  max_path_len_(max_path_len),
                  //TODO:: config
                  min_gap_(int(g.k() + 10)) {}

    virtual ~PathGapCloser() {}
};

class GapExtensionChooserFactory {
    const Graph &g_;
protected:
    const Graph& g() const {
        return g_;
    }
public:
    GapExtensionChooserFactory(const Graph &g): g_(g) {
    }

    virtual ~GapExtensionChooserFactory() {}
    virtual shared_ptr<ExtensionChooser> CreateChooser(const BidirectionalPath &original_path, size_t position) const = 0;
};

class ReadCloudGapExtensionChooserFactory : public GapExtensionChooserFactory {
    typedef shared_ptr <barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr;
    const ScaffoldingUniqueEdgeStorage unique_storage_;
    barcode_extractor_ptr extractor_;
public:
    ReadCloudGapExtensionChooserFactory(const Graph& g, const ScaffoldingUniqueEdgeStorage& unique_storage,
                                        barcode_extractor_ptr extractor) :
            GapExtensionChooserFactory(g), unique_storage_(unique_storage), extractor_(extractor) {}
    virtual ~ReadCloudGapExtensionChooserFactory() {}
    virtual shared_ptr<ExtensionChooser> CreateChooser(const BidirectionalPath& original_path, size_t position) const override {

        const EdgeId target_edge = FindUniqueAfterPosition(original_path, position);
//        const Gap original_gap = original_path.GapAt(position);
//        const size_t reserve = 1000;
//        auto dij = DijkstraHelper::CreateBoundedDijkstra(g, original_gap.gap + reserve);
//        dij.Run()
        auto extension_chooser_ptr = std::make_shared<ReadCloudGapExtensionChooser>(this->g(), extractor_, target_edge, unique_storage_);
        return extension_chooser_ptr;
    }

private:
    EdgeId FindUniqueAfterPosition(const BidirectionalPath& path, const size_t position) const {
        for (size_t i = position; i != path.Size(); ++i) {
            if (unique_storage_.IsUnique(path.At(i))) {
                return path.At(i);
            }
        }
        return EdgeId(0);
    }
};

class GapExtenderFactory {
public:
    virtual ~GapExtenderFactory() {}
    virtual shared_ptr<PathExtender> CreateExtender(const BidirectionalPath &original_path, size_t position) const = 0;
};

class SameExtenderFactory : public GapExtenderFactory {
    const shared_ptr<PathExtender> extender_;
public:
    SameExtenderFactory(shared_ptr<PathExtender> extender) : extender_(extender) {
    }

    shared_ptr<PathExtender> CreateExtender(const BidirectionalPath &, size_t) const override {
        return extender_;
    }
};

class SimpleExtenderFactory : public GapExtenderFactory {
    const conj_graph_pack &gp_;
    const GraphCoverageMap &cover_map_;
    const shared_ptr<GapExtensionChooserFactory> chooser_factory_;
    static const size_t MAGIC_LOOP_CONSTANT = 1000;
public:
    SimpleExtenderFactory(const conj_graph_pack &gp,
                          const GraphCoverageMap &cover_map,
                          const shared_ptr<GapExtensionChooserFactory> chooser_factory):
            gp_(gp),
            cover_map_(cover_map),
            chooser_factory_(chooser_factory) {
    }

    shared_ptr<PathExtender> CreateExtender(const BidirectionalPath &original_path, size_t position) const override {
        return make_shared<SimpleExtender>(gp_, cover_map_,
                                           chooser_factory_->CreateChooser(original_path, position),
                                           MAGIC_LOOP_CONSTANT,
                                           false,
                                           false);
    }
};

//Intermediate abstract class - majority of GapClosers needs only one next edge after gap, not entire original path.
class TargetEdgeGapCloser : public PathGapCloser {
protected:
    //returns updated gap to target edge
    virtual Gap CloseGap(EdgeId target_edge, const Gap &gap, BidirectionalPath &path) const = 0;

    Gap CloseGap(const BidirectionalPath &original_path,
                 size_t position, BidirectionalPath &path) const final override {
        return CloseGap(original_path.At(position), original_path.GapAt(position), path);
    }

public:
    TargetEdgeGapCloser(const Graph& g, size_t max_path_len):
            PathGapCloser(g, max_path_len) {}

};

class PathExtenderGapCloser: public TargetEdgeGapCloser {
    std::shared_ptr<path_extend::PathExtender> extender_;

protected:
    Gap CloseGap(const BidirectionalPath &original_path,
                 size_t position, BidirectionalPath &path) const override;

public:
    PathExtenderGapCloser(const Graph& g, size_t max_path_len, std::shared_ptr<PathExtender> extender):
            TargetEdgeGapCloser(g, max_path_len), extender_(extender) {
        DEBUG("ext added");
    }

    PathExtenderGapCloser(const Graph& g, size_t max_path_len,
                          shared_ptr<GapExtenderFactory> extender_factory):
            PathGapCloser(g, max_path_len),
            extender_factory_(extender_factory) {
        DEBUG("ext factory added");
    }

private:
    DECL_LOGGER("PathExtenderGapCloser")
};

class MatePairGapCloser: public TargetEdgeGapCloser {
    const std::shared_ptr<PairedInfoLibrary> lib_;
    const ScaffoldingUniqueEdgeStorage &storage_;
//TODO: config? somewhere else?
    static constexpr double weight_priority = 5;

    EdgeId FindNext(const BidirectionalPath &path,
                    const std::set<EdgeId> &present_in_paths,
                    VertexId last_v, EdgeId target_edge) const;
protected:
    Gap CloseGap(EdgeId target_edge, const Gap &gap, BidirectionalPath &path) const override;

    DECL_LOGGER("MatePairGapCloser")

public:
    MatePairGapCloser(const Graph &g, size_t max_path_len,
                      const std::shared_ptr<PairedInfoLibrary> lib,
                      const ScaffoldingUniqueEdgeStorage &storage):
            TargetEdgeGapCloser(g, max_path_len), lib_(lib), storage_(storage) {}
};

//TODO switch to a different Callback, no need to store all paths
class DijkstraGapCloser: public TargetEdgeGapCloser {
    typedef std::vector<std::vector<EdgeId>> PathsT;

    Gap FillWithMultiplePaths(const PathsT& paths,
                              BidirectionalPath& result) const;

    Gap FillWithBridge(const Gap &orig_gap,
                       const PathsT& paths, BidirectionalPath& result) const;

    size_t MinPathLength(const PathsT& paths) const;

    size_t MinPathSize(const PathsT& paths) const;

    std::vector<EdgeId> LCP(const PathsT& paths) const;

    std::map<EdgeId, size_t> CountEdgesQuantity(const PathsT& paths, size_t length_limit) const;

protected:
    Gap CloseGap(EdgeId target_edge, const Gap &gap, BidirectionalPath &path) const override;

    DECL_LOGGER("DijkstraGapCloser")

public:
    DijkstraGapCloser(const Graph& g, size_t max_path_len):
        TargetEdgeGapCloser(g, max_path_len) {}

};

class PathPolisher {
    static const size_t MAX_POLISH_ATTEMPTS = 5;

    const Graph &g_;
    std::vector<std::shared_ptr<PathGapCloser>> gap_closers_;

    void InfoAboutGaps(const PathContainer& result);

    std::unique_ptr<BidirectionalPath> Polish(const BidirectionalPath& path);
    DECL_LOGGER("PathPolisher")

public:
    PathPolisher(const Graph &g,
                 const std::vector<std::shared_ptr<PathGapCloser>> &gap_closers):
            g_(g), gap_closers_(gap_closers) {
    }

    PathContainer PolishPaths(const PathContainer &paths);
};


}
