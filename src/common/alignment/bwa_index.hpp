//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2016-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/mapping_path.hpp"

extern "C" {
struct bwaidx_s;
typedef struct bwaidx_s bwaidx_t;

struct mem_opt_s;
typedef struct mem_opt_s mem_opt_t;

struct mem_alnreg_s;
typedef struct mem_alnreg_s mem_alnreg_t;
typedef struct mem_alnreg_vs mem_alnreg_v;
};

namespace alignment {

class BWAIndex {
  public:
    enum class AlignmentMode {
        Default,
        Illumina,
        IntraCtg,
        PacBio,
        Ont2D,
        HiC
    };

    enum class RetainAlignments {
        Default,       // Default for the mode
        All,           // All, including secondary
        OnlyPrimary,   // Only primary
        QualityPrimary // Primary with mapq > 20
    };

    // bwaidx / memopt are incomplete below, therefore we need to outline ctor
    // and dtor.
    // FIXME: the retain logic should be bound to AlignSequence, but we're
    // making it state variable for the simplicity
    BWAIndex(const debruijn_graph::Graph& g,
             AlignmentMode mode = AlignmentMode::Default,
             RetainAlignments retain = RetainAlignments::Default);
    ~BWAIndex();

    omnigraph::MappingPath<debruijn_graph::EdgeId> AlignSequence(const Sequence &sequence,
                                                                 bool only_simple = false) const;
  private:
    void Init();
    omnigraph::MappingPath<debruijn_graph::EdgeId> GetMappingPath(const mem_alnreg_v&, const Sequence&, bool = false) const;

    const debruijn_graph::Graph& g_;

    // Store the options in memory
    std::unique_ptr<mem_opt_t, void(*)(void*)> memopt_;

    // hold the full index structure
    std::unique_ptr<bwaidx_t, void(*)(bwaidx_t*)> idx_;

    std::vector<debruijn_graph::EdgeId> ids_;

    AlignmentMode mode_;
    RetainAlignments retain_;

    DECL_LOGGER("BWAIndex");
};

}
