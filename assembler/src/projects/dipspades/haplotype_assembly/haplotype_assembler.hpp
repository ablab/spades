//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "conservative_regions_searcher.hpp"

namespace dipspades {

class HaplotypeAssembler {

    conj_graph_pack &consensus_graph_pack_;
    conj_graph_pack &double_graph_pack_;
    ContigStoragePtr default_storage_;
    ContigStoragePtr composite_storage_;
    CorrectionResult redundancy_map_;

public:

    HaplotypeAssembler(conj_graph_pack &consensus_graph_pack,
            conj_graph_pack &double_graph_pack,
            ContigStoragePtr default_storage,
            ContigStoragePtr composite_storage,
            CorrectionResult redundancy_map) :
                consensus_graph_pack_(consensus_graph_pack),
                double_graph_pack_(double_graph_pack),
                default_storage_(default_storage),
                composite_storage_(composite_storage),
                redundancy_map_(redundancy_map) {
            double_graph_pack_.kmer_mapper.Attach();
     }

    void Run() {
        INFO("Contigs separation starts");
        DiploidContigSeparator separator(consensus_graph_pack_.g, default_storage_,
                composite_storage_, redundancy_map_);
        INFO("Haplocontigs number: " << default_storage_->Size());
        INFO("Consensus contigs number: " << composite_storage_->Size());
        separator.SeparateContigs();
        SignedLabels signed_labels = separator.GetSignedLabels();
        string hapl_output(path::append_path(dsp_cfg::get().io.output_dir, "haplotype_assembly.out").c_str());
        signed_labels.WriteToFile(hapl_output, default_storage_);
        INFO("Result of haplotype assembly written in file " << hapl_output);
        INFO("Contigs separation ends");

        INFO("Conservative regions search starts");
        ConservativeRegionStorage conservative_regions = separator.GetConservativeRegionStorage();
        ConservativeRegionsSearcher cons_regions_searcher(double_graph_pack_, default_storage_,
                signed_labels, conservative_regions);
        cons_regions_searcher.Search();
        INFO("Conservative regions search ends");
    }
};

}
