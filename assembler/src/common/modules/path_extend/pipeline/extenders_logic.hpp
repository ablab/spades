//
// Created by andrey on 14.11.16.
//

#pragma once

#include "modules/path_extend/path_extender.hpp"
#include "launch_support.hpp"

namespace path_extend {

using namespace debruijn_graph;

struct ExtenderTriplet {
    io::LibraryType lib_type_;
    size_t lib_index_;
    shared_ptr<PathExtender> extender_;

    ExtenderTriplet(io::LibraryType lib_type, size_t lib_index, shared_ptr<PathExtender> extender):
        lib_type_(lib_type), lib_index_(lib_index), extender_(extender) {

    }

    bool operator<(const ExtenderTriplet& that) const {
        if (this->lib_type_ == that.lib_type_)
            return this->lib_index_ < that.lib_index_;
        return this->lib_type_ < that.lib_type_;
    }
};

typedef vector<ExtenderTriplet> ExtenderTriplets;

typedef vector<shared_ptr<PathExtender>> Extenders;

inline Extenders ExtractExtenders(const ExtenderTriplets& triplets) {
    Extenders result;
    for (const auto& triplet : triplets)
        result.push_back(triplet.extender_);

    return result;
}

class ExtendersGenerator {
    const config::dataset &dataset_info_;
    const PathExtendParamsContainer &params_;
    const conj_graph_pack &gp_;

    const GraphCoverageMap &cover_map_;
    const UniqueData &unique_data_;
    UsedUniqueStorage &used_unique_storage_;

    const PELaunchSupport &support_;

public:
    ExtendersGenerator(const config::dataset &dataset_info,
                       const PathExtendParamsContainer &params,
                       const conj_graph_pack &gp,
                       const GraphCoverageMap &cover_map,
                       const UniqueData &unique_data,
                       UsedUniqueStorage &used_unique_storage,
                       const PELaunchSupport& support) :
        dataset_info_(dataset_info),
        params_(params),
        gp_(gp),
        cover_map_(cover_map),
        unique_data_(unique_data),
        used_unique_storage_(used_unique_storage),
        support_(support) { }

    Extenders MakePBScaffoldingExtenders() const;

    Extenders MakeBasicExtenders() const;

    Extenders MakeMPExtenders() const;

    Extenders MakeCoverageExtenders() const;

    Extenders MakePEExtenders() const;

private:

    shared_ptr<SimpleExtender> MakePEExtender(size_t lib_index, bool investigate_loops) const;

    Extenders MakeMPExtenders(const ScaffoldingUniqueEdgeStorage &storage) const;

    shared_ptr<ExtensionChooser> MakeLongReadsExtensionChooser(size_t lib_index, const GraphCoverageMap& read_paths_cov_map) const;

    shared_ptr<SimpleExtender> MakeLongReadsExtender(size_t lib_index, const GraphCoverageMap& read_paths_cov_map) const;

    shared_ptr<SimpleExtender> MakeLongEdgePEExtender(size_t lib_index,
                                                      bool investigate_loops) const;

    shared_ptr<GapAnalyzer> MakeGapAnalyzer(double is_variation) const;

    shared_ptr<PathExtender> MakeScaffoldingExtender(size_t lib_index) const;

    shared_ptr<PathExtender> MakeRNAScaffoldingExtender(size_t lib_index) const;

    shared_ptr<PathExtender> MakeMatePairScaffoldingExtender
        (size_t lib_index, const ScaffoldingUniqueEdgeStorage &storage) const;

    shared_ptr<SimpleExtender> MakeCoordCoverageExtender(size_t lib_index) const;

    shared_ptr<SimpleExtender> MakeRNAExtender(size_t lib_index, bool investigate_loops) const;

    shared_ptr<SimpleExtender> MakeSimpleCoverageExtender(size_t lib_index) const;

    void PrintExtenders(const vector<shared_ptr<PathExtender>> &extenders) const;

};

}
