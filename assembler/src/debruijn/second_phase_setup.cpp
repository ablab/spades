//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "standard.hpp"
#include "dataset_readers.hpp"
#include "read_converter.hpp"

#include "de/paired_info.hpp"

#include "utils.hpp"
#include "stats/debruijn_stats.hpp"
#include "io/library.hpp"
#include "is_counter.hpp"
#include "pair_info_count.hpp"
#include "sequence_mapper.hpp"
#include "short_read_mapper.hpp"
#include "long_read_mapper.hpp"
#include "pair_info_filler.hpp"
#include "second_phase_setup.hpp"
#include "stats/debruijn_stats.hpp"
#include "path_extend/split_graph_pair_info.hpp"


namespace debruijn_graph {

void SecondPhaseSetup::run(conj_graph_pack &gp, const char*) {
	INFO("Preparing second phase");
	gp.ClearRRIndices();

	std::string old_pe_contigs_filename = cfg::get().output_dir + "scaffolds.fasta";
	std::string new_pe_contigs_filename = cfg::get().output_dir + "first_pe_contigs.fasta";

    VERIFY(path::check_existence(old_pe_contigs_filename));
    INFO("Moving preliminary contigs from " << old_pe_contigs_filename << " to " << new_pe_contigs_filename);
	int code = rename(old_pe_contigs_filename.c_str(), new_pe_contigs_filename.c_str());
    VERIFY(code == 0);

	io::ReadStreamList<io::SingleRead> additional_contigs;
	additional_contigs.push_back(io::EasyStream(cfg::get().output_dir + "first_pe_contigs.fasta", true));

	io::SequencingLibrary<debruijn_graph::debruijn_config::DataSetData> untrusted_contigs;
	untrusted_contigs.push_back_single(new_pe_contigs_filename);
	untrusted_contigs.set_orientation(io::LibraryOrientation::Undefined);
	untrusted_contigs.set_type(io::LibraryType::PathExtendContigs);
	cfg::get_writable().ds.reads.push_back(untrusted_contigs);

    //FIXME get rid of this awful variable
    cfg::get_writable().use_single_reads = false;
	INFO("Ready to run second phase");
}

}
