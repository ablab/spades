#pragma once

#include <memory>
#include <utility>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include "io/reads/paired_read.hpp"
#include "modules/assembly_graph/paths/mapping_path.hpp"
#include <modules/assembly_graph/graph_alignment/sequence_mapper.hpp>

using std::string;
using std::istringstream;

namespace tslr_resolver {
    typedef debruijn_graph::ConjugateDeBruijnGraph Graph;
    typedef runtime_k::RtSeq seq_t;
    typedef debruijn_graph::KmerFreeEdgeIndex<Graph, runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq>> KmerEdgeIndex;
    typedef debruijn_graph::EdgeIndex<Graph, seq_t, KmerEdgeIndex> Index;
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;
    typedef string BarcodeId;
    typedef std::unordered_set <BarcodeId> BarcodeSet;
    typedef std::unordered_map <EdgeId, BarcodeSet> barcode_map_t;
    typedef omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper;
    typedef debruijn_graph::KmerMapper<Graph> KmerSubs;
    typedef std::map <size_t, size_t> barcode_distribution;


    namespace tenx_barcode_parser {
        static const int barcode_len = 16;

        template <typename ReadType>
        bool is_valid(const ReadType &read) {
            auto str = read.name();
            return str.length() > barcode_len && str[barcode_len] == '#';
        }

        template <typename ReadType>
        BarcodeId GetTenxBarcode(const ReadType &read) {
            return read.name().substr(0, barcode_len);
        }
    } //tenx_barcode_parser

    struct tslr_barcode_library {
        string left_;
        string right_;
        string barcode_;
    };


    class BarcodeMapper {
    private:
        barcode_map_t barcode_map_heads;
        barcode_map_t barcode_map_tails;
        const Graph &g;
        const Index &index;
        const KmerSubs &kmer_mapper;
        size_t tail_threshold_;
        size_t norm_len = 10000;
    public:
        BarcodeMapper(const Graph &g, const Index& index,
                      const KmerSubs& kmer_mapper, size_t tail_threshold = 10000) :
                g(g), index(index), kmer_mapper(kmer_mapper), tail_threshold_(tail_threshold)
        {
            barcode_map_heads = barcode_map_t();
            barcode_map_tails = barcode_map_t();
        }

        BarcodeMapper(const BarcodeMapper& other) = default;

        void InitialFillMap(const Graph &g) {
            edge_it_helper helper(g);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                BarcodeSet set;
                barcode_map_heads.insert({*it, set});
                barcode_map_tails.insert({*it, set});
            }
        }

        void FillMap(const string &reads_filename, bool debug_mode = false) {
            auto lib_vec = GetLibrary(reads_filename);
            auto mapper = std::make_shared<debruijn_graph::NewExtendedSequenceMapper<Graph, Index> >
                          (g, index, kmer_mapper);

            int debug_counter = 0;

            for (auto lib: lib_vec) {
                std::string barcode = lib.barcode_;
                io::SeparatePairedReadStream paired_read_stream(lib.left_, lib.right_, 1);
                io::PairedRead read;
                while (!paired_read_stream.eof() && debug_counter < 100000) {
                    paired_read_stream >> read;
                    auto path_first = mapper -> MapRead(read.first());
                    auto path_second = mapper -> MapRead(read.second());
	            std::vector<omnigraph::MappingPath<EdgeId> > paths = {path_first, path_second};
		    for (auto path : paths) {	
                        for(size_t i = 0; i < path.size(); i++) {
                            if (is_at_edge_head(path[i].second))
                                barcode_map_heads.at(path[i].first).insert(barcode);
                            if (is_at_edge_tail(path[i].first, path[i].second))
                                barcode_map_tails.at(path[i].first).insert(barcode);
                        }
                    }
                    if (debug_mode) {
                        debug_counter++;
                    }
                }
            }
        }

        bool is_at_edge_tail(const EdgeId& edge, const omnigraph::MappingRange& range) {
            return range.mapped_range.start_pos + tail_threshold_ > g.length(edge);
        }

        bool is_at_edge_head(const omnigraph::MappingRange& range) {
            return range.mapped_range.end_pos < tail_threshold_;
        }

        BarcodeSet GetSetHeads(const EdgeId &edge) const {
            return barcode_map_heads.at(edge);
        }

        BarcodeSet GetSetTails(const EdgeId& edge) const {
            return barcode_map_tails.at(edge);
        }

        //Get barcode intersection of tail of first with head of second
        size_t IntersectionSize(const EdgeId &edge1, const EdgeId &edge2) const {
            size_t result = 0;
            auto Set1 = GetSetTails(edge1);
            auto Set2 = GetSetHeads(edge2);
            for (auto it = Set1.begin(); it != Set1.end(); ++it) {
                auto it2 = Set2.find(*it);
                if (it2 != Set2.end()) {
                    result++;
                }
            }
            return result;
        }

        double IntersectionSizeNormalized(const EdgeId &edge1, const EdgeId &edge2) const {
            return static_cast <double> (IntersectionSize(edge1, edge2)) / static_cast <double> (g.length(edge2) + norm_len);
        }

        //
        void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge, const std::string& which) {
            VERIFY(which == "head" || which == "tail");
            if (which == "head")
                barcode_map_heads[edge].insert(barcode);
            else
                barcode_map_tails[edge].insert(barcode);
        }

        barcode_map_t::const_iterator cbegin_heads() const noexcept {
            return barcode_map_heads.cbegin();
        }

        barcode_map_t::const_iterator cend_heads() const noexcept {
            return barcode_map_heads.cend();
        }

        barcode_map_t::const_iterator cbegin_tails() const noexcept {
            return barcode_map_tails.cbegin();
        }

        barcode_map_t::const_iterator cend_tails() const noexcept {
            return barcode_map_tails.cend();
        }

        size_t size(const std::string& which_end) const {
            VERIFY(which_end == "head" || which_end == "tail");
            if (which_end == "head")
                return barcode_map_heads.size();
            return barcode_map_tails.size();
        }

        std::pair <double, double> AverageBarcodeCoverage() {
            edge_it_helper helper(g);
            int64_t barcodes_overall_heads = 0;
            int64_t barcodes_overall_tails = 0;
            int64_t edges = 0;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                edges++;
                barcodes_overall_heads += barcode_map_heads.at(*it).size();
                barcodes_overall_tails += barcode_map_tails.at(*it).size();
            }
            DEBUG("heads: " + std::to_string(barcodes_overall_heads));
            DEBUG("tails: " + std::to_string(barcodes_overall_tails));
            DEBUG(edges);
            return make_pair(static_cast <double> (barcodes_overall_heads) / static_cast <double> (edges),
                             static_cast <double> (barcodes_overall_tails) / static_cast <double> (edges));
        }

//        barcode_distribution GetBarcodeDistribution(const EdgeId &edge)  {
//            barcode_distribution distr;
//            auto Set = GetSet(edge);
//            for (auto it = Set.begin(); it != Set.end(); ++it) {
//                distr[it -> size()]++;
//            }
//            return distr;
//        }
//
//        void SerializeBarcodeDistribution(std::ofstream& fout, const EdgeId& edge) {
//            fout << edge.int_id() << std::endl;
//            auto distr = GetBarcodeDistribution(edge);
//            for (auto it : distr) {
//                fout << it.first << ' ' << it.second << std::endl;
//            }
//            fout << std::endl;
//        }
//
//        void SerializeAllDistributions(std::ofstream& fout) {
//            edge_it_helper helper(g);
//            for(auto it = helper.begin(); it != helper.end(); ++it) {
//                SerializeBarcodeDistribution(fout, *it);
//            }
//        }
//
//        void SerializeOverallDistribution(std::ofstream& fout) {
//            barcode_distribution distr;
//            for (auto it = barcode_map_.begin(); it != barcode_map_.end(); ++it) {
//                for (auto elem)
//            }
//        }



        DECL_LOGGER("BarcodeMapper")

    private:
        std::vector <tslr_barcode_library> GetLibrary(const string& reads_filename) {
            std::vector <tslr_barcode_library> lib_vec;
            std::ifstream fin;
            fin.open(reads_filename);
            string line;
            while (getline(fin, line)) {
                if (!line.empty()) {
                    istringstream tmp_stream(line);
                    tslr_barcode_library lib;
                    tmp_stream >> lib.barcode_;
                    tmp_stream >> lib.left_;
                    tmp_stream >> lib.right_;
                    lib_vec.push_back(lib);
                }
            }
            return lib_vec;
        }
    };


} //tslr_resolver


