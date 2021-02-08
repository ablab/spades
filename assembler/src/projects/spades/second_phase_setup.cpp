//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "second_phase_setup.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/graph_iterators.hpp"
#include "utils/filesystem/path_helper.hpp"

#include <unordered_set>

namespace debruijn_graph {

void SecondPhaseSetup::run(GraphPack &gp, const char*) {
    INFO("Preparing second phase");
    gp.ClearRRIndices();
    gp.ClearPaths();

    std::string old_pe_contigs_filename = fs::append_path(cfg::get().output_dir, cfg::get().co.contigs_name + ".fasta");
    std::string new_pe_contigs_filename = fs::append_path(cfg::get().output_dir, "first_pe_contigs.fasta");

    VERIFY(fs::check_existence(old_pe_contigs_filename));
    INFO("Moving preliminary contigs from " << old_pe_contigs_filename << " to " << new_pe_contigs_filename);
    int code = rename(old_pe_contigs_filename.c_str(), new_pe_contigs_filename.c_str());
    VERIFY(code == 0);

    io::SequencingLibrary<config::LibraryData> untrusted_contigs;
    untrusted_contigs.push_back_single(new_pe_contigs_filename);
    untrusted_contigs.set_orientation(io::LibraryOrientation::Undefined);
    untrusted_contigs.set_type(io::LibraryType::PathExtendContigs);
    cfg::get_writable().ds.reads.push_back(untrusted_contigs);

    //FIXME get rid of this awful variable
    VERIFY(!cfg::get().use_single_reads);
    INFO("Ready to run second phase");
}

}
