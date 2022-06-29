//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "stats.hpp"
#include "reduce.hpp"

#include "common/io/graph/gfa_writer.hpp"
#include "common/io/graph/gfa_reader.hpp"
#include "common/toolchain/utils.hpp"

int main(int argc, char *argv[]) {
    toolchain::create_console_logger();

    CommandLineArguments command_line_arguments(argc, argv);

    gfa::GFAReader reader(command_line_arguments.GetGraphPath());
    debruijn_graph::Graph g(reader.k());
    reader.to_graph(g);

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
        gfa::GFAWriter writer(g, os);
        writer.WriteSegmentsAndLinks();

        INFO("Reduced graph is saved to " << output);
    }

    return 0;
}