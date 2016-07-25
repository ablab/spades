#pragma once

#include <memory>
#include <utility>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <bitset>
#include "io/reads/paired_read.hpp"
#include "modules/assembly_graph/paths/mapping_path.hpp"
#include <modules/assembly_graph/graph_alignment/sequence_mapper.hpp>

using std::string;
using std::istringstream;

namespace tslr_resolver {
    constexpr size_t max_barcodes = 384;


    typedef debruijn_graph::ConjugateDeBruijnGraph Graph;
    typedef runtime_k::RtSeq seq_t;
    typedef debruijn_graph::KmerFreeEdgeIndex<Graph, runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq>> KmerEdgeIndex;
    typedef debruijn_graph::EdgeIndex<Graph, seq_t, KmerEdgeIndex> Index;
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;
    typedef string BarcodeId;
    typedef std::bitset<max_barcodes> barcode_set_t;
    typedef std::unordered_map <EdgeId, barcode_set_t> barcode_map_t;
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

    class BarcodeEncoder {
        std::unordered_map <BarcodeId, size_t> codes_;
        size_t max_index;
    public:
        BarcodeEncoder() :
            codes_(), max_index (0)
        { }

        void AddEntry (const string& barcode) {
            auto it = codes_.find(barcode);
            if (it == codes_.end()) {
                codes_[barcode] = max_index;
            }
            VERIFY(max_index < max_barcodes);
            max_index++;
        }

        size_t GetCode (const string& barcode) const {
            VERIFY(codes_.find(barcode) != codes_.end());
            return codes_.at(barcode);
        }

        size_t GetSize() const {
            return max_index;
        }
    };

    class BarcodeMapper { 
        typedef string BarcodeId;
        //TODO: Make separate storing strategy
//        typedef std::unordered_set <BarcodeId> BarcodeSet;
    private:
        const Graph& g_;
        const Index& index_;
        const KmerSubs& kmer_mapper_;
        BarcodeEncoder barcode_codes_;
        barcode_map_t barcode_map_heads;
        barcode_map_t barcode_map_tails;
        size_t tail_threshold_;
        size_t norm_len_;
    public:
        BarcodeMapper(const Graph& g, const Index& index,
                      const KmerSubs& kmer_mapper, size_t tail_threshold = 10000, size_t norm_len = 10000) :
                g_(g), index_(index), kmer_mapper_(kmer_mapper), barcode_codes_(),
                tail_threshold_(tail_threshold), norm_len_(norm_len)
        {
            barcode_map_heads = barcode_map_t();
            barcode_map_tails = barcode_map_t();
        }

        void InitialFillMap(const Graph &g) {
            edge_it_helper helper(g);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                barcode_set_t set;
                barcode_map_heads.insert({*it, set});
                barcode_map_tails.insert({*it, set});
            }
        }

        void FillMap(const string &reads_filename, bool debug_mode = false) {
            auto lib_vec = GetLibrary(reads_filename);
            auto mapper = std::make_shared<debruijn_graph::NewExtendedSequenceMapper<Graph, Index> >
                          (g_, index_, kmer_mapper_);

            int debug_counter = 0;

            for (auto lib: lib_vec) {
                std::string barcode = lib.barcode_;
                barcode_codes_.AddEntry(barcode);
                size_t code = barcode_codes_.GetCode(barcode);
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
                                barcode_map_heads.at(path[i].first).set(code);
                            if (is_at_edge_tail(path[i].first, path[i].second))
                                barcode_map_tails.at(path[i].first).set(code);
                        }
                    }
                    if (debug_mode) {
                        debug_counter++;
                    }
                }
            }
        }

        bool is_at_edge_tail(const EdgeId& edge, const omnigraph::MappingRange& range) {
            return range.mapped_range.start_pos + tail_threshold_ > g_.length(edge);
        }

        bool is_at_edge_head(const omnigraph::MappingRange& range) {
            return range.mapped_range.end_pos < tail_threshold_;
        }

        barcode_set_t GetSetHeads(const EdgeId &edge) const {
            return barcode_map_heads.at(edge);
        }

        barcode_set_t GetSetTails(const EdgeId& edge) const {
            return barcode_map_tails.at(edge);
        }

        //Get barcode intersection of tail of first with head of second
        size_t IntersectionSize(const EdgeId &edge1, const EdgeId &edge2) const {
            auto Set1 = GetSetTails(edge1);
            auto Set2 = GetSetHeads(edge2);
            return (Set1 & Set2).count();
        }

        size_t UnionSize(const EdgeId &edge1, const EdgeId &edge2) const {
            auto Set1 = GetSetTails(edge1);
            auto Set2 = GetSetHeads(edge2);
            return (Set1 | Set2).count();
        }

        double IntersectionSizeNormalized(const EdgeId &edge1, const EdgeId &edge2) const {
            return static_cast <double> (IntersectionSize(edge1, edge2)) / static_cast <double> (g_.length(edge2) + norm_len_);
        }

        double IntersectionSizeRelative(const EdgeId &edge1, const EdgeId &edge2) const {
            return static_cast <double> (IntersectionSize(edge1, edge2)) / static_cast <double> (UnionSize(edge1, edge2));
        }

        //
        void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge, const std::string& which) {
            VERIFY(which == "head" || which == "tail");
            size_t code = barcode_codes_.GetCode(barcode);
            if (which == "head")
                barcode_map_heads[edge].set(code);
            else
                barcode_map_tails[edge].set(code);
        }

        void InsertSet (barcode_set_t& set, const EdgeId& edge, const std::string& which) {
            VERIFY(which == "head" || which == "tail");
            VERIFY(set.size() == max_barcodes);
            if (which == "head")
                barcode_map_heads[edge] = set;
            else
                barcode_map_tails[edge] = set;
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

        size_t size() const {
            return barcode_map_tails.size();
        }

        std::pair <double, double> AverageBarcodeCoverage() {
            edge_it_helper helper(g_);
            int64_t barcodes_overall_heads = 0;
            int64_t barcodes_overall_tails = 0;
            int64_t edges = 0;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                edges++;
                barcodes_overall_heads += barcode_map_heads.at(*it).count();
                barcodes_overall_tails += barcode_map_tails.at(*it).count();
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


