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
#include <modules/pipeline/config_struct.hpp>

using std::string;
using std::istringstream;

namespace tslr_resolver {
    constexpr int16_t max_barcodes = 384;

    typedef debruijn_graph::ConjugateDeBruijnGraph Graph;
    typedef runtime_k::RtSeq seq_t;
    typedef debruijn_graph::KmerFreeEdgeIndex<Graph, runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq>> KmerEdgeIndex;
    typedef debruijn_graph::EdgeIndex<Graph, seq_t, KmerEdgeIndex> Index;
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;
    typedef omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper;
    typedef debruijn_graph::KmerMapper<Graph> KmerSubs;
    typedef string BarcodeId;


    struct tslr_barcode_library {
        string left_;
        string right_;
        string barcode_;
    };

    class BarcodeEncoder {
        std::unordered_map <BarcodeId, int16_t> codes_;
        int16_t max_index;
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

        int16_t GetCode (const string& barcode) const {
            VERIFY(codes_.find(barcode) != codes_.end());
            return codes_.at(barcode);
        }

        int16_t GetSize() const {
            return max_index;
        }
    };

    class BarcodeMapper {
    public:
    protected:
        const Graph& g_;
    public:
        BarcodeMapper (const Graph &g) :
                g_(g) {}
        virtual ~BarcodeMapper() {}
        virtual size_t size() const = 0;
        virtual void FillMap (const string& reads_filename, const Index& index,
                              const KmerSubs& kmer_mapper) = 0;
        virtual size_t IntersectionSize(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual size_t UnionSize(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual double IntersectionSizeNormalizedByUnion(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual double IntersectionSizeNormalizedBySecond(const EdgeId &edge1, const EdgeId &edge2) const = 0;
        virtual double IntersectionSizeNormalizedByFirst(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual double AverageBarcodeCoverage () const = 0;
        virtual size_t GetSizeHeads(const EdgeId& edge) const = 0;
        virtual size_t GetSizeTails(const EdgeId& edge) const = 0;
        virtual void ReadEntry(ifstream& fin, const EdgeId& edge) = 0;
        virtual void WriteEntry(ofstream& fin, const EdgeId& edge) = 0;
        virtual void FilterByAbundance(size_t threshold) = 0;
        virtual void SerializeOverallDistribution(const string& path) const = 0;

    };

    template <class barcode_entry_t>
    class HeadTailBarcodeMapper : public BarcodeMapper {
    protected:
        typedef std::unordered_map <EdgeId, barcode_entry_t> barcode_map_t;
        using BarcodeMapper::g_;
        barcode_map_t barcode_map_heads;
        size_t tail_threshold_;
        BarcodeEncoder barcode_codes_;

    public:
        HeadTailBarcodeMapper (const Graph &g, size_t tail_threshold) :
                BarcodeMapper(g), barcode_map_heads(),
                tail_threshold_(tail_threshold),  barcode_codes_() {
            InitialFillMap();
        }

        HeadTailBarcodeMapper (const HeadTailBarcodeMapper& other) = default;

        virtual ~HeadTailBarcodeMapper() {}

        void InitialFillMap() {
            edge_it_helper helper(g_);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                barcode_entry_t set(*it);
                barcode_map_heads.insert({*it, set});
            }
        }

        size_t size() const {
            return barcode_map_heads.size();
        }

        typename barcode_map_t::const_iterator cbegin() const noexcept {
            return barcode_map_heads.cbegin();
        }

        typename barcode_map_t::const_iterator cend() const noexcept {
            return barcode_map_heads.cend();
        }


        virtual double IntersectionSizeNormalizedByUnion(const EdgeId& edge1, const EdgeId& edge2) const override {
            return static_cast <double> (IntersectionSize(edge1, edge2)) / 
                static_cast <double> (UnionSize(edge1, edge2));
        }

        virtual double IntersectionSizeNormalizedBySecond(const EdgeId &edge1, const EdgeId &edge2) const override {
            return static_cast <double> (IntersectionSize(edge1, edge2)) / 
                static_cast <double> (GetSizeHeads(edge2));
        }

        virtual double IntersectionSizeNormalizedByFirst(const EdgeId &edge1, const EdgeId& edge2) const override {
            return static_cast <double> (IntersectionSize(edge1, edge2)) / 
                static_cast <double> (GetSizeTails(edge1));
        }


        size_t GetSizeHeads(const EdgeId& edge) const override {
            return barcode_map_heads.at(edge).Size();
        }

        size_t GetSizeTails(const EdgeId& edge) const override {
            return barcode_map_heads.at(g_.conjugate(edge)).Size();
        }

        void FillMap (const string& reads_filename, const Index& index,
                      const KmerSubs& kmer_mapper) {
            auto lib_vec = GetLibrary(reads_filename);
            auto mapper = std::make_shared<debruijn_graph::NewExtendedSequenceMapper<Graph, Index> >
                    (g_, index, kmer_mapper);

            for (auto lib: lib_vec) {
                std::string barcode = lib.barcode_;
                barcode_codes_.AddEntry(barcode);
                std::shared_ptr<io::ReadStream<io::PairedRead>> paired_read_ptr =
                        make_shared<io::SeparatePairedReadStream> (lib.left_, lib.right_, 1);
                auto wrapped_stream = io::RCWrap(paired_read_ptr);
                io::PairedRead read;
                while (!wrapped_stream->eof()) {
                    *wrapped_stream >> read;
                    auto path_first = mapper -> MapRead(read.first());
                    auto path_second = mapper -> MapRead(read.second());
                    std::vector<omnigraph::MappingPath<EdgeId> > paths = {path_first, path_second};
                    for (auto path : paths) {
                        for(size_t i = 0; i < path.size(); i++) {
                            if (is_at_edge_head(path[i].second))
                                InsertBarcode(barcode, path[i].first);
                            if (is_at_edge_tail(path[i].first, path[i].second))
                                InsertBarcode(barcode, g_.conjugate(path[i].first));
                        }
                    }
                }
            }
        }

        double AverageBarcodeCoverage() const override {
            edge_it_helper helper(g_);
            int64_t barcodes_overall_heads = 0;
            int64_t barcodes_overall_tails = 0;
            int64_t long_edges = 0;
            size_t len_threshold = cfg::get().ts_res.len_threshold;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                if (g_.length(*it) > len_threshold) {
                    long_edges++;
                    barcodes_overall_heads += GetSizeHeads(*it);
                    barcodes_overall_tails += GetSizeTails(*it);
                }
            }
            DEBUG("heads: " + std::to_string(barcodes_overall_heads));
            DEBUG("tails: " + std::to_string(barcodes_overall_tails));
            DEBUG("Long edges" + long_edges);
            return static_cast <double> (barcodes_overall_heads) / static_cast <double> (long_edges);
        }

        size_t IntersectionSize(const EdgeId &edge1, const EdgeId &edge2) const override {
            return barcode_map_heads.at(edge1).IntersectionSize(edge2);
        }

        size_t UnionSize(const EdgeId &edge1, const EdgeId &edge2) const override {
            return barcode_map_heads.at(edge1).UnionSize(edge2);
        }


        //Delete low abundant barcodes from every edge
        void FilterByAbundance(size_t trimming_threshold) override {
            for (auto entry = barcode_map_heads.begin(); entry != barcode_map_heads.end(); ++entry) {
                entry->second.Filter(trimming_threshold);
            }
        }

        void SerializeOverallDistribution(const string& path) const override {
            ofstream fout;
            fout.open(path);
            std::map <size_t, size_t> overall_distr;
            for (auto entry: barcode_map_heads) {
                auto current_distr = barcode_map_heads.at(entry.first);
                for (auto it = current_distr.GetDistribution().cbegin();
                     it != current_distr.GetDistribution().cend() ; ++it) {
                    overall_distr[it->second]++;
                }
            }
            for (auto entry : overall_distr) {
                fout << entry.first << ": " << entry.second << endl;
            }
        }

        void ReadEntry (ifstream& fin, const EdgeId& edge) override {
            barcode_entry_t entry(edge);
            size_t distr_size;
            fin >> distr_size;
            for (size_t i = 0; i < distr_size; ++i) {
                int16_t bid;
                size_t abundance;
                fin >> bid >> abundance;
                entry.InsertBarcode(bid, abundance);
            };
            barcode_map_heads[edge] = entry;
            DEBUG(edge.int_id());
            DEBUG(distr_size);
        }

        void WriteEntry (ofstream& fout, const EdgeId& edge) override {
            fout << g_.int_id(edge) << std::endl;
            auto distribution = GetEntryHeads(edge).GetDistribution();
            fout << distribution.size() << endl;
            for (auto entry : distribution) {
                fout << entry.first << ' ' << entry.second << endl;
            }
        }

    protected:
        barcode_entry_t GetEntryHeads(const EdgeId &edge) {
            return barcode_map_heads.at(edge);
        }

        barcode_entry_t GetEntryTails(const EdgeId &edge) {
            return barcode_map_heads.at(g_.conjugate(edge));
        }

        void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge) {
            int16_t code = barcode_codes_.GetCode(barcode);
            barcode_map_heads.at(edge).InsertBarcode(code);
        }

        //utils
        bool is_at_edge_tail(const EdgeId& edge, const omnigraph::MappingRange& range) {
            return range.mapped_range.start_pos + tail_threshold_ > g_.length(edge);
        }

        bool is_at_edge_head(const omnigraph::MappingRange& range) {
            return range.mapped_range.end_pos < tail_threshold_;
        }

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


    class SimpleBarcodeEntry {
        typedef std::unordered_map <int64_t, size_t> barcode_distribution_t;
        EdgeId edge_;
        barcode_distribution_t barcode_distribution_;

    public:
        SimpleBarcodeEntry():
            edge_(), barcode_distribution_() {};
        SimpleBarcodeEntry(const EdgeId& edge) :
                edge_(edge), barcode_distribution_() {}

        SimpleBarcodeEntry(const SimpleBarcodeEntry& other) = default;
        SimpleBarcodeEntry& operator =(const    SimpleBarcodeEntry& other) = default;

        barcode_distribution_t GetDistribution() const {
            return barcode_distribution_;
        }

        size_t IntersectionSize(const SimpleBarcodeEntry& other) const {
            size_t result = 0;
            auto distr_this = barcode_distribution_;
            auto distr_other = other.GetDistribution();
            for (auto it = distr_this.begin(); it != distr_this.end(); ++it) {
                if (distr_other.find(it-> first) != distr_other.end()) {
                    result++;
                }
            }
            return result;
        }

        size_t UnionSize(const SimpleBarcodeEntry& other) const {
            auto distr_this = barcode_distribution_;
            auto distr_other = other.GetDistribution();
            return Size() + other.Size() - IntersectionSize(other);
        }

        void InsertBarcode(int16_t code, size_t abundance = 1) {
            if (barcode_distribution_.find(code) == barcode_distribution_.end()) {
                barcode_distribution_.insert({code, abundance});
            }
            else {
                barcode_distribution_.at(code) += abundance;
            }
        }

        void InsertSet (barcode_distribution_t& set) {
            barcode_distribution_ = set;
        }

        void Filter (size_t trimming_threshold) {
            for (auto it = barcode_distribution_.begin(); it != barcode_distribution_.end() ;) {
                if (it->second < trimming_threshold) {
                    DEBUG("Erased " + it->first + ' ' + std::to_string(it->second));
                    barcode_distribution_.erase(it++);
                }
                else {
                    ++it;
                }
            }
        }

        size_t Size() const {
            return barcode_distribution_.size();
        }

    };



} //tslr_resolver


