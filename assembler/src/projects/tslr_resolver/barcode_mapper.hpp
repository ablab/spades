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

    struct barcode_library {
        string left_;
        string right_;
        string barcode_;
    };


    class BarcodeMapper {
    private:
        barcode_map_t barcode_map_;
        const Graph &g;
        const Index &index;
        const KmerSubs &kmer_mapper;
    public:
        BarcodeMapper(const Graph &g, const Index& index,
                      const KmerSubs& kmer_mapper) :
                g(g), index(index), kmer_mapper(kmer_mapper)
        {
            barcode_map_ = barcode_map_t();
        }

        BarcodeMapper(const BarcodeMapper& other) = default;

        void FillMap(const string &reads_filename, bool debug_mode = false) {
                        edge_it_helper helper(g);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                BarcodeSet set;
                barcode_map_.insert({*it, set});
            }
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
                    for(size_t i = 0; i < path_first.size(); i++) {
                        barcode_map_.at(path_first[i].first).insert(barcode);
                    }
                    for(size_t i = 0; i < path_second.size(); i++) {
                        barcode_map_.at(path_second[i].first).insert(barcode);
                    }
                    if (debug_mode) {
                        debug_counter++;
                    }
                }
            }
        }

        BarcodeSet GetSet(const EdgeId &edge) const {
            return barcode_map_.at(edge);
        }

        size_t IntersectionSize(const EdgeId &edge1, const EdgeId &edge2) const {
            size_t result = 0;
            auto Set1 = GetSet(edge1);
            auto Set2 = GetSet(edge2);
            for (auto it = Set1.begin(); it != Set1.end(); ++it) {
                auto it2 = Set2.find(*it);
                if (it2 != Set2.end()) {
                    result++;
                }
            }
            return result;
        }

        void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge) {
            barcode_map_[edge].insert(barcode);
        }

        barcode_map_t::const_iterator cbegin() const noexcept {
            return barcode_map_.cbegin();
        }
        barcode_map_t::const_iterator cend() const noexcept {
            return barcode_map_.cend();
        }

        size_t size() const {
            return barcode_map_.size();
        }

        double AverageBarcodeCoverage() {
            edge_it_helper helper(g);
            int64_t barcodes_overall = 0;
            int64_t edges = 0;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                edges++;
                barcodes_overall += barcode_map_.at(*it).size();
            }
            INFO(barcodes_overall);
            INFO(edges);
            return static_cast <double> (barcodes_overall) / static_cast <double> (edges);
        }

        barcode_distribution GetBarcodeDistribution(const EdgeId &edge)  {
            barcode_distribution distr;
            auto Set = GetSet(edge);
            for (auto it = GetSet(edge).begin(); it != GetSet(edge).end(); ++it) {
                distr[it -> size()]++;
            }
            return distr;
        }

        void SerializeBarcodeDistribution(std::ofstream& fout, const EdgeId& edge) {
            fout << edge.int_id() << std::endl;
            for (auto it : GetBarcodeDistribution(edge)) {
                fout << it.first << ' ' << it.second << std::endl;
            }
        }

        DECL_LOGGER("BarcodeMapper")

    private:
        std::vector <barcode_library> GetLibrary(const string& reads_filename) {
            std::vector <barcode_library> lib_vec;
            std::ifstream fin;
            fin.open(reads_filename);
            string line;
            while (getline(fin, line)) {
                if (!line.empty()) {
                    istringstream tmp_stream(line);
                    barcode_library lib;
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


