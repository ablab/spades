//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "find_scope.hpp"
#include "stats.hpp"
#include "reduce.hpp"

#include "common/io/graph/gfa_writer.hpp"
#include "common/io/graph/gfa_reader.hpp"
#include "common/toolchain/utils.hpp"

int main(int argc, char *argv[]) {
    toolchain::create_console_logger();

    CommandLineArguments command_line_arguments(argc, argv);

    io::IdMapper<std::string> *id_mapper = new io::IdMapper<std::string>();
    gfa::GFAReader reader(command_line_arguments.GetGraphPath());
    debruijn_graph::Graph g(reader.k());
    reader.to_graph(g, id_mapper);

    if (command_line_arguments.GetCommand() == "stat") {
        gfa_tools::CalculateStat(g, command_line_arguments.GetNPercentiles(),
                                 command_line_arguments.GetLengthPercentile(),
                                 command_line_arguments.GetSequencesCoveragePercentiles(),
                                 command_line_arguments.GetYamlOutputPath());
    }

    if (command_line_arguments.GetCommand() == "reduce") {
        gfa_tools::MakeReduction(g, command_line_arguments);

        std::filesystem::path output = command_line_arguments.GetGraphOutputPath();
        std::ofstream os(output);
        gfa::GFAWriter writer(g, os, io::MapNamingF<debruijn_graph::Graph>(*id_mapper));
        writer.WriteSegmentsAndLinks();

        INFO("Reduced graph is saved to " << output);
    }

    if (command_line_arguments.GetCommand() == "around_sequence_scope") {
            std::string sequence_name = command_line_arguments.GetSequenceNameToFindScope();
            debruijn_graph::EdgeId edge_id((*id_mapper)[sequence_name]);
            size_t scope_depth = command_line_arguments.GetScopeDepth();
            omnigraph::GraphComponent<debruijn_graph::Graph> scope = 
                                gfa_tools::FindAroundSequenceScope(g, edge_id, scope_depth);
            std::filesystem::path output = command_line_arguments.GetGraphOutputPath();
            std::ofstream os(output);
            gfa::GFAComponentWriter writer(scope, os, io::MapNamingF<debruijn_graph::Graph>(*id_mapper));
            writer.WriteSegmentsAndLinks();
            INFO("Scope around sequence " << sequence_name << " is saved to " << output);
    }
    return 0;  
}