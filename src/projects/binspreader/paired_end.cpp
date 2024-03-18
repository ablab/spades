//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "paired_end.hpp"

#include "link_index.hpp"

#include "assembly_graph/core/graph.hpp"

#include "paired_info/paired_info_utils.hpp"
#include "alignment/kmer_sequence_mapper.hpp"

#include "io/binary/paired_index.hpp"
#include "io/dataset_support/read_converter.hpp"

#include <threadpool/threadpool.hpp>

using namespace debruijn_graph;

namespace binning {

void FillPairedEndLinks(LinkIndex &pe_links,
                        SequencingLib &lib,
                        const debruijn_graph::Graph &graph,
                        const std::filesystem::path &workdir,
                        unsigned nthreads,
                        bool bin_load, bool bin_save) {
    if (!bin_load || !io::ReadConverter::LoadLibIfExists(lib)) {
        std::unique_ptr<ThreadPool::ThreadPool> pool;
        if (nthreads > 1)
            pool = std::make_unique<ThreadPool::ThreadPool>(nthreads);
        io::ReadConverter::ConvertToBinary(lib, pool.get());
    }

    paired_info::PairedIndex index(graph);
    if (!bin_load) {
        alignment::ShortKMerReadMapper mapper(graph, workdir);

        paired_info::FillPairedIndex(graph,
                                     mapper,
                                     lib, index, { }, 0, std::numeric_limits<unsigned>::max());

        if (bin_save) {
            INFO("Saving paired-end information");
            io::binary::Save(workdir / "paired_index", index);
        }
    } else {
        INFO("Loading paired-end information");
        io::binary::Load(workdir / "paired_index", index);
    }
              
    for (EdgeId e1 : graph.edges()) {
        for (auto entry : index.GetHalf(e1)) {
            EdgeId e2 = entry.first, ce2 = graph.conjugate(e2);
            VERIFY(entry.second.size() == 1);
            if (!(e2 <= ce2))
                continue;

            // Exclude self-links
            if (e1 == e2)
                continue;
                      
            omnigraph::de::DEWeight w = 0;
            auto AddToWeight = [&](EdgeId e1, EdgeId e2) {
                const auto& hist = index.Get(e1, e2);
                if (!hist.empty())
                    w += hist.begin()->weight;
            };

            if (e2 != e1) {
                // Need to aggregate links:
                // e1 => e2, e1 => e2', e2' => e1, e2 => e1
                AddToWeight(e1, e2); AddToWeight(e1, ce2);
                AddToWeight(e2, e1); AddToWeight(ce2, e1);
            } else { // e2 == e1
                // Need to aggregate links:
                // e1 => e1, e1' => e1, e1 => e1'
                AddToWeight(e1, e1); AddToWeight(ce2, e1); AddToWeight(e1, ce2);
            }

            pe_links.add(e1, e2, w);
        }
    }
}

}
