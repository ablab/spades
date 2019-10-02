//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/stats/picture_dump.hpp"
#include <io/reads/osequencestream.hpp>
#include <common/io/reads/binary_converter.hpp>
#include <common/io/reads/edge_sequences_reader.hpp>
#include "assembly_graph/components/connected_component.hpp"
#include "assembly_graph/stats/statistics.hpp"
#include "assembly_graph/paths/path_finders.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include <common/utils/filesystem/path_helper.hpp>

namespace debruijn_graph {

inline void OutputEdgeSequences(const Graph &g, const std::string &contigs_output_filename) {
    INFO("Outputting contigs to " << contigs_output_filename << ".fasta");
    io::osequencestream_cov oss(contigs_output_filename + ".fasta");

    for (EdgeId e: g.canonical_edges()) {
        oss << g.coverage(e);
        oss << g.EdgeNucls(e).str();
    }
}

inline void OutputEdgeSequencesToBinary(const Graph &g, const std::string &contigs_output_dir) {
    INFO("Outputting contigs to " << contigs_output_dir);
    fs::make_dir(contigs_output_dir);

    io::BinaryWriter single_converter(contigs_output_dir + "/contigs");
    io::SingleStream single_reader = io::EdgeSequencesStream(g);
    single_converter.ToBinary(single_reader);
}

inline void OutputEdgesByID(const Graph &g,
                                const std::string &contigs_output_filename) {
    INFO("Outputting contigs to " << contigs_output_filename << ".fasta");
    io::OFastaReadStream oss(contigs_output_filename + ".fasta");
    for (EdgeId e: g.canonical_edges()) {
        std::string s = g.EdgeNucls(e).str();
        oss << io::SingleRead(io::MakeContigId(g.int_id(e), s.size(), g.coverage(e), "EDGE"), s);
    }
}
} // namespace debruijn_grap

