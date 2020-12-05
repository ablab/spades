//***************************************************************************
//* Copyright (c) 2016-2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "modules/path_extend/path_extender.hpp"
#include "modules/path_extend/gap_analyzer.hpp"
#include "launch_support.hpp"

namespace path_extend {

using namespace debruijn_graph;

struct ExtenderTriplet {
    io::LibraryType lib_type_;
    size_t lib_index_;
    std::shared_ptr<PathExtender> extender_;

    ExtenderTriplet(io::LibraryType lib_type, size_t lib_index, std::shared_ptr<PathExtender> extender):
        lib_type_(lib_type), lib_index_(lib_index), extender_(extender) {

    }

    static int GetPriority(io::LibraryType type) {
        #define SET_PRIORITY(o) case o: --priority
        using io::LibraryType;
        int priority = 0;
        // the higher the element, the lower the return value
        switch (type) {
            // would have the minimum return value
            SET_PRIORITY(LibraryType::SangerReads);
            SET_PRIORITY(LibraryType::PacBioReads);
            SET_PRIORITY(LibraryType::NanoporeReads);
            SET_PRIORITY(LibraryType::TrustedContigs);
            SET_PRIORITY(LibraryType::SingleReads);
            SET_PRIORITY(LibraryType::PairedEnd);
            SET_PRIORITY(LibraryType::HQMatePairs);
            SET_PRIORITY(LibraryType::MatePairs);
            SET_PRIORITY(LibraryType::TSLReads);
            SET_PRIORITY(LibraryType::PathExtendContigs);
            SET_PRIORITY(LibraryType::UntrustedContigs);
            SET_PRIORITY(LibraryType::FLRNAReads);
            SET_PRIORITY(LibraryType::AssemblyGraph);
            // would have the maximum return value
            break;
            // there is a LibraryType that does not have priority
            default: VERIFY(false);
        };
        return priority;
        #undef SET_PRIORITY
    }

    bool operator<(const ExtenderTriplet& that) const {
        if (GetPriority(this->lib_type_) == GetPriority(that.lib_type_))
            return this->lib_index_ < that.lib_index_;
        return GetPriority(this->lib_type_) < GetPriority(that.lib_type_);
    }
};

typedef std::vector<ExtenderTriplet> ExtenderTriplets;

typedef std::vector<std::shared_ptr<PathExtender>> Extenders;

inline Extenders ExtractExtenders(const ExtenderTriplets& triplets) {
    Extenders result;
    for (const auto& triplet : triplets)
        result.push_back(triplet.extender_);

    return result;
}

class ExtendersGenerator {
    const config::dataset &dataset_info_;
    const PathExtendParamsContainer &params_;
    const GraphPack &gp_;
    const Graph &graph_;

    const GraphCoverageMap &cover_map_;
    const UniqueData &unique_data_;
    UsedUniqueStorage &used_unique_storage_;

    const PELaunchSupport &support_;

public:
    ExtendersGenerator(const config::dataset &dataset_info,
                       const PathExtendParamsContainer &params,
                       const GraphPack &gp,
                       const GraphCoverageMap &cover_map,
                       const UniqueData &unique_data,
                       UsedUniqueStorage &used_unique_storage,
                       const PELaunchSupport& support) :
        dataset_info_(dataset_info),
        params_(params),
        gp_(gp),
        graph_(gp.get<Graph>()),
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

    std::shared_ptr<SimpleExtender> MakePEExtender(size_t lib_index, bool investigate_loops) const;

    Extenders MakeMPExtenders(const ScaffoldingUniqueEdgeStorage &storage) const;

    std::shared_ptr<ExtensionChooser> MakeLongReadsExtensionChooser(size_t lib_index,
                                                                    const GraphCoverageMap &read_paths_cov_map) const;

    std::shared_ptr<SimpleExtender> MakeLongReadsExtender(size_t lib_index,
                                                          const GraphCoverageMap &read_paths_cov_map) const;

    std::shared_ptr<SimpleExtender> MakeLongEdgePEExtender(size_t lib_index,
                                                      bool investigate_loops) const;

    std::shared_ptr<GapAnalyzer> MakeGapAnalyzer(double is_variation) const;

    std::shared_ptr<PathExtender> MakeScaffoldingExtender(size_t lib_index) const;

    std::shared_ptr<PathExtender> MakeRNAScaffoldingExtender(size_t lib_index) const;

    std::shared_ptr<PathExtender> MakeMatePairScaffoldingExtender(size_t lib_index,
                                                                  const ScaffoldingUniqueEdgeStorage &storage) const;

    std::shared_ptr<SimpleExtender> MakeCoordCoverageExtender(size_t lib_index) const;

    std::shared_ptr<SimpleExtender> MakeRNAExtender(size_t lib_index, bool investigate_loops) const;

    std::shared_ptr<SimpleExtender> MakeSimpleCoverageExtender(size_t lib_index) const;

    std::shared_ptr<ExtensionChooser> MakeLongReadsRNAExtensionChooser(size_t lib_index, const GraphCoverageMap& read_paths_cov_map) const;

    std::shared_ptr<SimpleExtender> MakeLongReadsRNAExtender(size_t lib_index, const GraphCoverageMap& read_paths_cov_map) const;

    void PrintExtenders(const std::vector<std::shared_ptr<PathExtender>> &extenders) const;

};

}
