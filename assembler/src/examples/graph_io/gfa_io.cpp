//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gfa_io.hpp"

#include "common/io/graph/gfa_reader.hpp"
#include "common/io/graph/gfa_writer.hpp"

/*
 * GFA (Graphical Fragment Assembly) is a file format for fragment assembly. It represents deBruijn graph
 * as a set of segments (edges) and links between segments (vertexes). Each segment has a unique label. File also
 * stores information about segments.
 */

namespace spades_example {

void SaveToGFA(const debruijn_graph::Graph& g, const std::filesystem::path& save_to="saved_graph.gfa") {
    /*
     * GFAWriter is declared in common/io/graph/gfa_writer.hpp file. To utilize it user should include this
     * file and link graphio library to the project.
     *
     * An instance of GFAWriter class should be initialised by 2 fields: a reference to the graph which will be
     * saved in GFA format, and an ostream instance, in which graph will be written. In this example we write
     * graph to a file, hence we transmit an ofstream instance.
     */
    std::ofstream os(save_to);
    gfa::GFAWriter writer(g, os);
    /*
     * The WriteSegmentsAndLinks() function is writing a graph in GFA format into the output stream.
     */
    writer.WriteSegmentsAndLinks();
}

void ReadFromGFA(debruijn_graph::Graph& g, const std::filesystem::path& read_from) {
    /*
     * GFAReader is declared in common/io/graph/gfa_reader.hpp file. To utilize it user should include this
     * file and link graphio library to the project.
     *
     * An instance of GFAReader class should be initialised by a path to a file in GFA format.
     */
    gfa::GFAReader reader(read_from);
    /*
     * The to_graph(debruijn_graph::Graph& g) function takes a non-const reference to a debruijn_graph::Graph
     * class instance. The graph from GFA file will be written in this instance.
     */
    reader.to_graph(g);
}
} // spades_example