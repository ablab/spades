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

    enum MapperType {
        Bitset,
        Trimmable
    };

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
        MapperType type_;
    protected:
        const Graph& g_;
    public:
        BarcodeMapper (const Graph &g, MapperType type) :
                type_(type), g_(g) {}
        virtual ~BarcodeMapper() {}
        virtual size_t size() const = 0;
        virtual void FillMap (const string& reads_filename, const Index& index,
                              const KmerSubs& kmer_mapper, bool debug_mode = false) = 0;
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
    public:
        using BarcodeMapper::type_;
    protected:
        typedef std::unordered_map <EdgeId, barcode_entry_t> barcode_map_t;
        using BarcodeMapper::g_;
        barcode_map_t barcode_map_heads;
        size_t tail_threshold_;
        size_t norm_len_;
        BarcodeEncoder barcode_codes_;

    public:
        HeadTailBarcodeMapper (const Graph &g, MapperType type, size_t tail_threshold, size_t norm_len) :
                BarcodeMapper(g, type), barcode_map_heads(),
                tail_threshold_(tail_threshold), norm_len_(norm_len), barcode_codes_() {
            InitialFillMap();
        }

        HeadTailBarcodeMapper (const HeadTailBarcodeMapper& other) = default;

        virtual ~HeadTailBarcodeMapper() {}

        void InitialFillMap() {
            edge_it_helper helper(g_);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                barcode_entry_t set;
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
                static_cast <double> (GetSizeHeads(edge2) * GetSizeHeads(edge2));
        }

        virtual double IntersectionSizeNormalizedByFirst(const EdgeId &edge1, const EdgeId& edge2) const override {
            return static_cast <double> (IntersectionSize(edge1, edge2)) / 
                static_cast <double> (GetSizeTails(edge1));
        }

        barcode_entry_t GetSetHeads(const EdgeId &edge) const {
            return barcode_map_heads.at(edge);
        }

        barcode_entry_t GetSetTails(const EdgeId& edge) const {
            return barcode_map_heads.at(g_.conjugate(edge));
        }

        size_t GetSizeHeads(const EdgeId& edge) const override {
            return GetSetHeads(edge).size();
        }

        size_t GetSizeTails(const EdgeId& edge) const override {
            return GetSetTails(edge).size();
        }

        void FillMap (const string& reads_filename, const Index& index,
                      const KmerSubs& kmer_mapper, bool debug_mode = false) {
            auto lib_vec = GetLibrary(reads_filename);
            auto mapper = std::make_shared<debruijn_graph::NewExtendedSequenceMapper<Graph, Index> >
                    (g_, index, kmer_mapper);

            int debug_counter = 0;
            int debug_steps = 100000;

            for (auto lib: lib_vec) {
                std::string barcode = lib.barcode_;
                barcode_codes_.AddEntry(barcode);
                io::SeparatePairedReadStream paired_read_stream(lib.left_, lib.right_, 1);
                io::PairedRead read;
                while (!paired_read_stream.eof() && debug_counter < debug_steps) {
                    paired_read_stream >> read;
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
                    if (debug_mode) {
                        debug_counter++;
                        if (debug_counter % 10000 == 0) {
                            INFO(std::to_string(debug_counter) + " reads processed...");
                        }
                    }
                }
            }
        }

    protected:
        barcode_entry_t GetSetHeads(const EdgeId &edge) {
            return barcode_map_heads.at(edge);
        }

        barcode_entry_t GetSetTails(const EdgeId& edge) {
            return barcode_map_heads.at(g_.conjugate(edge));
        }

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

        virtual void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge) = 0;

        virtual void InsertSet (barcode_entry_t& set, const EdgeId& edge) = 0;
    };


    class BitSetBarcodeMapper : public HeadTailBarcodeMapper<std::bitset<max_barcodes> > {
        typedef string BarcodeId;
        typedef std::bitset<max_barcodes> barcode_set_t;
    public:
        using HeadTailBarcodeMapper::type_;
    private:
        using HeadTailBarcodeMapper::g_;
        using HeadTailBarcodeMapper::barcode_map_heads;
        using HeadTailBarcodeMapper::barcode_codes_;
    public:
        BitSetBarcodeMapper(const Graph& g, MapperType type,
                            size_t tail_threshold = 10000, size_t norm_len = 10000) :
                HeadTailBarcodeMapper(g, type, tail_threshold, norm_len)
        {}

        size_t IntersectionSize(const EdgeId &edge1, const EdgeId &edge2) const override {
            auto Set1 = GetSetTails(edge1);
            auto Set2 = GetSetHeads(edge2);
            return (Set1 & Set2).count();
        }

        size_t UnionSize(const EdgeId &edge1, const EdgeId &edge2) const override {
            auto Set1 = GetSetTails(edge1);
            auto Set2 = GetSetHeads(edge2);
            return (Set1 | Set2).count();
        }

        double AverageBarcodeCoverage() const override {
            edge_it_helper helper(g_);
            int64_t barcodes_overall_heads = 0;
            int64_t barcodes_overall_tails = 0;
            int64_t edges = 0;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                edges++;
                barcodes_overall_heads += barcode_map_heads.at(*it).count();
            }
            DEBUG("heads: " + std::to_string(barcodes_overall_heads));
            DEBUG("tails: " + std::to_string(barcodes_overall_tails));
            DEBUG(edges);
            return static_cast <double> (barcodes_overall_heads) / static_cast <double> (edges);
        }

        void ReadEntry (ifstream& fin, const EdgeId& edge) override {
            barcode_set_t entry;
            fin >> entry;
            InsertSet(entry, edge);
        }

        void WriteEntry (ofstream& fout, const EdgeId& edge) override {
            barcode_set_t set;
            set = GetSetHeads(edge);
            fout << edge << ' ' << set << std::endl;
        }

        //FIXME Stupid way to avoid warning. Filtering will be moved to FillMap to allow compressing
        void FilterByAbundance(size_t trimming_threshold) override {
            VERIFY(trimming_threshold > 1);
        }

        void SerializeOverallDistribution(const string& path) const override {
            ofstream fout;
            fout.open(path);
        }

    private:

        void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge) override {
            int16_t code = barcode_codes_.GetCode(barcode);
            barcode_map_heads[edge].set(code);
        }

        void InsertSet (barcode_set_t& set, const EdgeId& edge) override {
            barcode_map_heads[edge] = set;
        }
        DECL_LOGGER("BarcodeMapperBitset")
    };


    class TrimmableBarcodeMapper : public HeadTailBarcodeMapper<std::unordered_map <int16_t, size_t> > {
        typedef string BarcodeId;
        typedef std::unordered_map <int16_t, size_t> barcode_distribution_t;
    public:
        using HeadTailBarcodeMapper::type_;
    private:
        using HeadTailBarcodeMapper::g_;
        using HeadTailBarcodeMapper::barcode_map_heads;
        using HeadTailBarcodeMapper::barcode_codes_;
    public:
        TrimmableBarcodeMapper(const Graph& g, MapperType type,
                               size_t tail_threshold = 10000, size_t norm_len = 10000) :
                HeadTailBarcodeMapper(g, type, tail_threshold, norm_len)
        {}

        size_t IntersectionSize(const EdgeId &edge1, const EdgeId &edge2) const override {
            size_t result = 0;
            auto Distr1 = GetSetTails(edge1);
            auto Distr2 = GetSetHeads(edge2);
            for (auto it = Distr1.begin(); it != Distr1.end(); ++it) {
                if (Distr2.find(it-> first) != Distr2.end()) {
                    result++;
                }
            }
            return result;
        }

        size_t UnionSize(const EdgeId &edge1, const EdgeId &edge2) const override {
            auto Distr1 = GetSetTails(edge1);
            auto Distr2 = GetSetHeads(edge2);
            return Distr1.size() + Distr2.size() - IntersectionSize(edge1, edge2);
        }


        double AverageBarcodeCoverage() const override {
            edge_it_helper helper(g_);
            int64_t barcodes_overall_heads = 0;
            int64_t barcodes_overall_tails = 0;
            int64_t edges = 0;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                edges++;
                barcodes_overall_heads += barcode_map_heads.at(*it).size();
                barcodes_overall_tails += barcode_map_heads.at(g_.conjugate(*it)).size();
            }
            DEBUG("heads: " + std::to_string(barcodes_overall_heads));
            DEBUG("tails: " + std::to_string(barcodes_overall_tails));
            DEBUG(edges);
            return static_cast <double> (barcodes_overall_heads) / static_cast <double> (edges);
        }

        //Delete low abundant barcodes
        void FilterByAbundance(size_t trimming_threshold) override {
            for (auto entry = barcode_map_heads.begin(); entry != barcode_map_heads.end(); ++entry) {
                for (auto it = entry->second.begin(); it != entry->second.end() ;) {
                    if (it->second < trimming_threshold) {
                        DEBUG("Erased " + it->first + ' ' + std::to_string(it->second));
                        entry->second.erase(it++);
                    }
                    else {
                        ++it;
                    }
                }
            }
        }

        void SerializeOverallDistribution(const string& path) const override {
            ofstream fout;
            fout.open(path);
            std::map <size_t, size_t> overall_distr;
            for (auto entry: barcode_map_heads) {
                const barcode_distribution_t& current_distr = barcode_map_heads.at(entry.first);
                for (auto it = current_distr.cbegin(); it != current_distr.cend() ; ++it) {
                    overall_distr[it->second]++;
                }
            }
            for (auto entry : overall_distr) {
                fout << entry.first << ": " << entry.second << endl;
            }
        }

        void ReadEntry (ifstream& fin, const EdgeId& edge) override {
            barcode_distribution_t distr;
            size_t distr_size;
            fin >> distr_size;
            for (size_t i = 0; i < distr_size; ++i) {
                int16_t bid;
                size_t abundance;
                fin >> bid >> abundance;
                distr.insert({bid, abundance});
            };
            barcode_map_heads[edge] = distr;
            DEBUG(edge.int_id());
            DEBUG(distr_size);
        }

        void WriteEntry (ofstream& fout, const EdgeId& edge) override {
            barcode_distribution_t distr;
            fout << g_.int_id(edge) << std::endl;
            distr = GetSetHeads(edge);
            fout << distr.size() << endl;
            for (auto entry : distr) {
                fout << entry.first << ' ' << entry.second << endl;
            }
        }

    private:
        void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge) override {
            int16_t code = barcode_codes_.GetCode(barcode);
            if (barcode_map_heads.at(edge).find(code) == barcode_map_heads.at(edge).end()) {
                barcode_map_heads.at(edge).insert({code, 1});
            }
            else {
                barcode_map_heads.at(edge).at(code)++;
            }
        }

        void InsertSet (barcode_distribution_t& set, const EdgeId& edge) override {
                barcode_map_heads[edge] = set;
        }

        DECL_LOGGER("TrimmableBarcodeMapper")
    };


} //tslr_resolver


