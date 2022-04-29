//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "common/assembly_graph/core/construction_helper.hpp"
#include "common/assembly_graph/core/graph.hpp"
#include "common/utils/logger/log_writers.hpp"
#include "examples/graph_io/gfa_io.hpp"

#include <random>
#include <unordered_map>

#include <clipp/clipp.h>

/*
 * debruijn_graph::EdgeId - the unique identifier of an edge in the graph.
 * debruijn_graph::VertexId - the unique identifier of a vertex in a graph.
 *
 * The information about the edge is stored in debruijn_graph::DeBruijnEdgeData class
 * (common/assembly_graph/core/debruijn_data.hpp). Every edge has the following properties: sequence
 * of nucleotides (class Sequence), coverage (debruijn_graph::CoverageData),
 * and flanking coverage (debruijn_graph::CoverageData).
 *
 * Inside the debruijn_graph::CoverageData coverage stores in uint64_t.
 *
 * Class Sequence stores the sequence of nucleotides. Roughly, it is a string of chars, where nucleotides
 * are represented as 'A' for adenine, 'C' for cytosine, 'G' for guanine, and 'T' for thymine. An instance of
 * Sequence can be constructed from string of ACGT / acgt / 0123 chars.
 *
 *
 * To get the sequence of an edge of the graph g with EdgeId = e, use g.EdgeNucls(e) method.
 * To get std::string from Sequence s, use s.str() method.
 * To get an average coverage (double) of an edge of the graph g with EdgeId = e, use g.coverage(e).
 * All these methods are declared in assembly_graph/core/graph.cpp.
 *
 * Vertexes of debruijn_graph::Graph do not store information about nucleotide sequence, but it is possible
 * to get it via g.VertexNucls(v) method, where g is a debruijn_graph::Graph, and v is a
 * debruijn_graph::VertexId.
 *
 */


void AddBamboo(debruijn_graph::Graph& g) {
    INFO("Start of adding new edges and vertexes to create a bamboo in a graph")
   /*
    * This is an artificial example of adding extra edges and vertexes to a graph to create a bamboo to
    * use the ReplaceBamboo(debruijn_graph::Graph& g) function. New edges here are not correct according
    * to rules of deBruijn graph construction.
    *
    * Class omnigraph::ConstructionHelper allows to clone the graph, add new edges, create new vertexes,
    * and link edges with vertexes and edges with each other. Class is declared in
    * common/assembly_graph/core/construction_helper.hpp.
    *
    * To get a ConstructionHelper of the graph g, use g.GetConstructionHelper() method.
    */

    omnigraph::ConstructionHelper helper = g.GetConstructionHelper();

    auto it = g.SmartVertexBegin();
    debruijn_graph::VertexId a = *it;
    ++it;
    debruijn_graph::VertexId b = *(it);

    /*
     * To add a new edge via ConstructionHelper user should create DeBruijnEdgeData instance first. The
     * constructor of this class takes a sequence of nucleotides (class Sequence).
     *
     * ConstructionHelper::AddEdge(DeBruijnEdgeData edata) method returns an EdgeId of a new edge.
     */
    debruijn_graph::DeBruijnEdgeData edata1(Sequence("ATTGAACC"));
    debruijn_graph::EdgeId e1 = helper.AddEdge(edata1);

    debruijn_graph::DeBruijnEdgeData edata2(Sequence("TTGAACCC"));
    debruijn_graph::EdgeId e2 = helper.AddEdge(edata2);

    debruijn_graph::DeBruijnEdgeData edata3(Sequence("TGAACCCC"));
    debruijn_graph::EdgeId e3 = helper.AddEdge(edata3);

    /*
     * To add a new vertex to a graph use ConstructionHelper::AddVertex(DeBruijnVertexData vdata) method.
     * It can take a new DeBruijnVertexData() rvalue. This method returns a VertexId of a new vertex.
     */
    debruijn_graph::VertexId v1 = helper.CreateVertex(debruijn_graph::DeBruijnVertexData());
    debruijn_graph::VertexId v2 = helper.CreateVertex(debruijn_graph::DeBruijnVertexData());

    /*
     * To link a vertex with outgoing or incoming edge use ConstructionHelper::LinkOutgoingEdge(VertexId, EdgeId)
     * or ConstructionHelper::LinkIncomingEdge(VertexId, EdgeId) methods, respectively.
     */
    helper.LinkOutgoingEdge(a, e1);
    helper.LinkIncomingEdge(v1, e1);

    helper.LinkOutgoingEdge(v1, e2);
    helper.LinkIncomingEdge(v2, e2);

    helper.LinkOutgoingEdge(v2, e3);
    helper.LinkIncomingEdge(b, e3);

    INFO("Bamboo edges are added into the graph");
}


void ReplaceBamboo(debruijn_graph::Graph& g) {
    INFO("Start of the bamboo replacement");

    omnigraph::ConstructionHelper helper = g.GetConstructionHelper(); // assembly_graph/core/construction_helper.hpp

    /*
     * There are several iterators to go through vertexes in debruijn_graph::Graph. They are declared in
     * common/assembly_graph/core/observable_graph.hpp file.
     *
     * Vertex iterators return a pointer to debruijn_graph::VertexId.
     */
    for (auto it = g.SmartVertexBegin(); !it.IsEnd(); ++it) {
        /*
         * Methods g.OutgoingEdgeCount(v) and g.IncomingEdgeCount(v) return number of outgoing and
         * incoming edges of the vertex with debruijn_graph::VertexId v respectively. These methods are declared
         * in common/assembly_graph/core/graph_core.hpp file.
         */
        if (g.OutgoingEdgeCount(*it) == 1 && g.IncomingEdgeCount(*it) == 1) {
            INFO("Vertex " << *it << " has one incoming and one outgoing edges");

            /*
             * Methods g.in_begin(v) and d.out_begin(v) return pointers to a first edges in a list of
             * incoming and outgoing edges of a vertex with VertexId v, respectively. In this case there is
             * only one incoming and outgoing edges, hence we can be sure that these pointers point
             * to the edges we need.
             * These methods are declared in common/assembly_graph/core/graph_core.hpp file.
             */
            debruijn_graph::EdgeId e1 = *g.in_begin(*it);
            debruijn_graph::EdgeId e2 = *g.out_begin(*it);

            /*
             * To get a VertexId of a vertex from or in which an edge with EdgeId e of graph g comes, use
             * g.EdgeStart(e) or g.EdgeEnd(e), respectively.
             */
            debruijn_graph::VertexId v1 = g.EdgeStart(e1);
            debruijn_graph::VertexId v2 = g.EdgeEnd(e2);

            /*
             * New edge should contain an intersection of nucleotide sequences from the replaced edges.
             * The length of overlap equals number of nucleotides in a vertex, which is k.
             */
            size_t k = g.k();
            //k = 7; // uncomment it, if you are replacing a bamboo, added via AddBamboo(g) function
            debruijn_graph::DeBruijnEdgeData edata(Sequence(g.EdgeNucls(e1) + g.EdgeNucls(e2).Subseq(k)));
            debruijn_graph::EdgeId new_edge = helper.AddEdge(edata);
            INFO("A new edge id created");

            /*
             * Method g.DeleteEdge(e) deletes edge with debruijn_graph::EdgeId e from graph g.
             * Method g.DeleteVertex(v) deletes vertex with debruijn_graph::VertexId from graph g.
             * These methods are declared in assembly_graph/core/observable_graph.hpp.
             *
             * It is important to delete edges linked to a vertex (or link them to other vertexes)
             * before deleting the vertex.
             */
            g.DeleteEdge(e1);
            g.DeleteEdge(e2);
            g.DeleteVertex(*it);

            INFO("Edges and vertexes of a bamboo are deleted");

            helper.LinkOutgoingEdge(v1, new_edge);
            helper.LinkIncomingEdge(v2, new_edge);

            INFO("A replacing edge is linked to vertexes");
        }
    }

    //uncomment it to check if all bamboos were replaced by single edges
    /*
    for (auto it = g.SmartVertexBegin(); !it.IsEnd(); ++it) {
        if (g.OutgoingEdgeCount(*it) == 1 && g.IncomingEdgeCount(*it) == 1) {
            INFO("Vertex " << *it << " has one incoming and one outgoing edges");
        }
    }
    */
    INFO("Bamboo edges are replaced with single edges");
}

// TODO: clip tips, remove one of the bubble arcs

void DeleteEdgesWithSomeFeature(debruijn_graph::Graph& g, double coverage_threshold) {
    /*
     * This function demonstrates the deletion of edges with some particular features from the graph.
     * In a case of this example edges with coverage lower than threshold and with the number of
     * some nucleotide higher, then the half of the edge length are deleted.

     * There are several iterators to go through edges in debruijn_graph::Graph. They are declared in
     * common/assembly_graph/core/observable_graph.hpp file.
     *
     * Edge iterators return a pointer to debruijn_graph::EdgeId.
     */
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        std::string sequence_ = g.EdgeNucls(*it).str();
        uint64_t nucl_cnt_ = std::max({std::count(sequence_.begin(), sequence_.end(), 'A'),
                                       std::count(sequence_.begin(), sequence_.end(), 'T'),
                                       std::count(sequence_.begin(), sequence_.end(), 'G'),
                                       std::count(sequence_.begin(), sequence_.end(), 'C')});
        double average_coverage_ = g.coverage(*it);
        if (nucl_cnt_ > sequence_.size() / 2 || average_coverage_ < coverage_threshold) {
            /*
             * Method g.DeleteEdge(e) deletes edge with debruijn_graph::EdgeId e from the graph g.
             * This method is declared in assembly_graph/core/observable_graph.hpp.
             */
             g.DeleteEdge(*it);
        }
    }
    // uncomment it to check if all edges were filtered
    /*
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
        std::cout << *it << " " << g.EdgeNucls(*it).str() << " " << g.coverage(*it) << "\n";
    }
    */
}

void CreateConsoleLogger(const std::filesystem::path& log_fn="") {
    using namespace logging;
    logger *lg = create_logger(exists(log_fn) ? log_fn : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

void ParseCommandLine(int argc, char *argv[], std::filesystem::path& load_graph_from, double& coverage_threshold) {
    using namespace clipp;
    std::string load_from_;
    auto cli = (
            option("-l", "--load_graph_from") & value("file", load_from_),
            option("-c", "--coverage_threshold") & opt_value("threshold", coverage_threshold)
    );
    parse(argc, argv, cli);

    load_graph_from = load_from_;
}


int main(int argc, char *argv[]) {
    CreateConsoleLogger();
    INFO("Start of the graph transformations example");

    std::filesystem::path load_graph_from;
    double coverage_threshold = 1;
    ParseCommandLine(argc, argv, load_graph_from, coverage_threshold);

    debruijn_graph::Graph g(55);
    spades_example::ReadFromGFA(g, load_graph_from);
    INFO("Graph is loaded");

    AddBamboo(g); // FIXME: how to call bamboo properly?
    ReplaceBamboo(g);
    DeleteEdgesWithSomeFeature(g, coverage_threshold);

    INFO("Graph transformation example is finished");
    return 0;
}