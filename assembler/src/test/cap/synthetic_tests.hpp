#pragma once

#include "compare_standard.hpp"
#include "logger/log_writers.hpp"
#include "graphio.hpp"
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
    typedef graph_pack<GraphT, Seq, DeBruijnEdgeIndex<KmerStoringDeBruijnEdgeIndex<GraphT, Seq>>> GraphPackT;

    const string filename_;
    const size_t k_;
    const string output_dir_;
    const string work_dir_;

    XmlTreeT xml_tree_;

    void ProcessExample(ContigStreamsPtr streams, size_t id) const {
        GraphPackT gp(k_, work_dir_);
        ColorHandler<GraphT> coloring(gp.g);
        ConstructColoredGraph(gp, coloring, *RCWrapStreams(*streams), /*fill_pos*/true);
        Save(gp, coloring, output_dir_ + ToString(id));
    }

    void Save(const GraphPackT& gp, const ColorHandler<GraphT>& coloring, const string& filename) const {
        typename PrinterTraits<GraphT>::Printer printer(gp.g, gp.int_ids);
        INFO("Saving graph to " << filename);
        printer.saveGraph(filename);
        printer.saveEdgeSequences(filename);
        printer.savePositions(filename, gp.edge_pos);
        SaveColoring(gp.g, gp.int_ids, coloring, filename);
        PrintColoredGraphWithColorFilter(gp.g, coloring, gp.edge_pos,
                                         filename + ".dot");
    }

    const vector<io::SingleRead> ParseGenome(
            const XmlTreeT& genome_node) const {
        vector<io::SingleRead> contigs;
        size_t contig_cnt = 0;
        BOOST_FOREACH(const XmlNodeT& contig, genome_node) {
            contigs.push_back(
                    io::SingleRead("contig_" + ToString(contig_cnt++),
                                   contig.second.data()));
        }
        return contigs;
    }

    void ProcessExample(const XmlTreeT& example_node, size_t id) const {

        vector<vector<io::SingleRead>> genomes;
        BOOST_FOREACH (const XmlNodeT& genome, example_node) {
            if (genome.first == "genome") {
                genomes.push_back(ParseGenome(genome.second));
            }
        }

        INFO("--------------------------------------------");
        INFO("Processing example " << id);
        VERIFY(genomes.size() == 2);
        ContigStreamsPtr streams(new ContigStreams());
        streams->push_back(new io::VectorReader<Contig>(genomes[0]));
        streams->push_back(new io::VectorReader<Contig>(genomes[1]));
        ProcessExample(streams, id);
    }

    vector<size_t> ProcessFile(const set<size_t>& ids_to_launch) const {
        vector<size_t> launched_tests;
        BOOST_FOREACH (const XmlNodeT& example, xml_tree_.get_child("examples")) {
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

}
;

}
