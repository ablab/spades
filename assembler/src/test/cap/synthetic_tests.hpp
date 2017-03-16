//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "compare_standard.hpp"
#include "comparison_utils.hpp"
#include "utils/logger/log_writers.hpp"
#include "pipeline/graphio.hpp"
#include <boost/test/unit_test.hpp>

#include "comparison_utils.hpp"
#include "diff_masking.hpp"
#include "repeat_masking.hpp"
#include "assembly_compare.hpp"
#include "test_utils.hpp"

namespace cap {

template<class Seq>
class SyntheticTestsRunner {
    typedef boost::property_tree::ptree XmlTreeT;
    typedef XmlTreeT::value_type XmlNodeT;
    typedef ConjugateDeBruijnGraph GraphT;
    typedef graph_pack<GraphT, Seq, KmerStoringEdgeIndex<GraphT, Seq, kmer_index_traits<Seq>, SimpleStoring>> GraphPackT;

    const string filename_;
    const size_t k_;
    const string output_dir_;
    const string work_dir_;

    XmlTreeT xml_tree_;

    size_t NumberOfContigs(ContigStream& stream) const {
        size_t cnt = 0;
        while(!stream.eof()) {
            Contig c;
            stream >> c;
            cnt++;
        }
        stream.reset();
        return cnt;
    }

    vector<string> TransparentContigNames(ContigStreams& streams) const {
        vector<string> genome_names;
        for (auto &stream : streams) {
            stream.reset();

            io::SingleRead contig;
            while (!stream.eof()) {
                stream >> contig;
                genome_names.push_back(contig.name());
            }

            stream.reset();
        }
        return genome_names;
    }

    void PrintBlocks(BlockPrinter<GraphT>& block_printer, ContigStreams& streams) const {
        vector<string> contig_names = TransparentContigNames(streams);

        size_t transparent_id = 0;
        for (size_t i = 0; i < streams.size(); ++i) {
            io::RCRemovingWrapper<Contig> stream(streams.ptr_at(i));
            for (size_t j = 0,
                    n = NumberOfContigs(stream);
                    j < n; ++j) {
                block_printer.ProcessContig(unsigned(i + 1), unsigned(transparent_id), contig_names[transparent_id]);
                transparent_id += 2;
            }
        }
    }

    void ProcessExample(ContigStreams streams, size_t id) const {
        GraphPackT gp(k_, work_dir_, 0);
        ColorHandler<GraphT> coloring(gp.g);
        CoordinatesHandler<GraphT> coordinates_handler;

        ConstructColoredGraph(gp, coloring, coordinates_handler, streams);
        coordinates_handler.FindGenomeFirstEdge(3);
        Save(gp, coloring, coordinates_handler, streams, output_dir_ + std::to_string(id));
    }

    void Save(const GraphPackT& gp, const ColorHandler<GraphT>& coloring,
            const CoordinatesHandler<GraphT> &coordinates_handler,
            ContigStreams& streams, const string& file_name) const {
        typename debruijn_graph::graphio::PrinterTraits<GraphT>::Printer printer(gp.g);
        INFO("Saving graph to " << file_name);
        printer.SaveGraph(file_name);
        printer.SaveEdgeSequences(file_name);
        //        printer.savePositions(filename, gp.edge_pos);
        SaveColoring(gp.g, coloring, file_name);

        shared_ptr<GraphSplitter<Graph>> splitter = omnigraph::ReliableSplitter(gp.g,
                numeric_limits<size_t>::max(),
                numeric_limits<size_t>::max());

//        LengthIdGraphLabeler<Graph> /*length_*/labeler(gp.g);
        LengthGraphLabeler<Graph> /*length_*/labeler(gp.g);
//        EdgeCoordinatesGraphLabeler<Graph> pos_labeler(gp.g,
//                                                       coordinates_handler,
//                                                       TransparentContigNames(*streams));
//
//        CompositeLabeler<Graph> labeler(length_labeler, pos_labeler);

        visualization::visualization_utils::WriteComponents(gp.g, file_name, splitter,
                coloring.ConstructColorer(), labeler);

        BlockPrinter<GraphT> block_printer(gp.g, coordinates_handler, file_name + ".blk");
        PrintBlocks(block_printer, streams);
        streams.reset();
    }

    const vector<io::SingleRead> ParseGenome(
            const XmlTreeT& genome_node) const {
        vector<io::SingleRead> contigs;
        size_t contig_cnt = 0;
        for (const XmlNodeT& contig : genome_node) {
            contigs.push_back(
                    io::SingleRead("contig_" + std::to_string(contig_cnt++),
                            contig.second.data()));
        }
        return contigs;
    }

    void ProcessExample(const XmlTreeT& example_node, size_t id) const {

        vector<vector<io::SingleRead>> genomes;
        for (const XmlNodeT& genome : example_node) {
            if (genome.first == "genome") {
                genomes.push_back(ParseGenome(genome.second));
            }
        }

        INFO("--------------------------------------------");
        INFO("Processing example " << id);
        VERIFY(genomes.size() == 2);
        ContigStreams streams;
        streams.push_back(make_shared<io::VectorReadStream<Contig>>(genomes[0]));
        streams.push_back(make_shared<io::VectorReadStream<Contig>>(genomes[1]));
        ProcessExample(io::RCWrap(streams), id);
    }

    vector<size_t> ProcessFile(const set<size_t>& ids_to_launch) const {
        vector<size_t> launched_tests;
        for (const XmlNodeT& example : xml_tree_.get_child("examples")) {
            size_t id = example.second.get<size_t>("<xmlattr>.n");
            if (ids_to_launch.empty() || ids_to_launch.count(id) > 0) {
                ProcessExample(example.second, id);
                launched_tests.push_back(id);
            }
        }
        return launched_tests;
    }

public:
    SyntheticTestsRunner(const string& filename, size_t k,
            const string& output_dir, const string& work_dir)
    : filename_(filename),
    k_(k),
    output_dir_(output_dir),
    work_dir_(work_dir) {
        read_xml(filename, xml_tree_);
    }

    vector<size_t> Run(const vector<size_t>& test_numbers = vector<size_t>()) {
        return ProcessFile(set<size_t>(test_numbers.begin(), test_numbers.end()));
    }

};

}
