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

    template <class barcode_entry_t>
    class BarcodeMapper { //TODO: Make separate HeadTailBarcodeMapper
    protected:
        typedef std::unordered_map <EdgeId, barcode_entry_t> barcode_map_t;
        const Graph &g_;
        barcode_map_t barcode_map_heads;
        barcode_map_t barcode_map_tails;

    public:
        MapperType type_;

    public:
        BarcodeMapper (const Graph &g, MapperType type) :
                g_(g), type_(type), barcode_map_heads(), barcode_map_tails() {
            InitialFillMap(g);
        }

        BarcodeMapper (const BarcodeMapper& other) = default;

        virtual ~BarcodeMapper() {}

        void InitialFillMap(const Graph &g) {
            edge_it_helper helper(g);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                barcode_entry_t set;
                barcode_map_heads.insert({*it, set});
                barcode_map_tails.insert({*it, set});
            }
        }

        size_t size() const {
            return barcode_map_heads.size();
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

        virtual void FillMap (const string& reads_filename, const Index& index,
                              const KmerSubs& kmer_mapper, bool debug_mode = false) = 0;
        //Get barcode intersection of tail of first with head of second
        virtual size_t IntersectionSize(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        //Get barcode union of tail of first with head of second
        virtual size_t UnionSize(const EdgeId& edge1, const EdgeId& edge2) const = 0;
        virtual std::pair <double, double> AverageBarcodeCoverage () const = 0;
        virtual void ReadEntry(ifstream& fin, const EdgeId& edge, const string& which_end) = 0;
        virtual void WriteEntry(ofstream& fin, const EdgeId& edge, const string& which_end) = 0;

        virtual double IntersectionSizeRelative(const EdgeId& edge1, const EdgeId& edge2) {
            return static_cast <double> (IntersectionSize(edge1, edge2)) / static_cast <double> (UnionSize(edge1, edge2));
        }
    };


    class BitSetBarcodeMapper : public BarcodeMapper<std::bitset<max_barcodes> > {
        typedef string BarcodeId;
        typedef std::bitset<max_barcodes> barcode_set_t;
    private:
        using BarcodeMapper::g_;
        using BarcodeMapper::type_;
        using BarcodeMapper::barcode_map_heads
        using BarcodeMapper::barcode_map_tails;
        BarcodeEncoder barcode_codes_;
        size_t tail_threshold_;
        size_t norm_len_;
    public:
        BitSetBarcodeMapper(const Graph& g, MapperType type,
                            size_t tail_threshold = 10000, size_t norm_len = 10000) :
                BarcodeMapper(g, type), barcode_codes_(),
                tail_threshold_(tail_threshold), norm_len_(norm_len)
        {
            barcode_map_heads = barcode_map_t();
            barcode_map_tails = barcode_map_t();
        }



        void FillMap(const string &reads_filename, const Index& index,
                     const KmerSubs& kmer_mapper, bool debug_mode = false) override {
            auto lib_vec = GetLibrary(reads_filename);
            auto mapper = std::make_shared<debruijn_graph::NewExtendedSequenceMapper<Graph, Index> >
                          (g_, index, kmer_mapper);

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

        double IntersectionSizeNormalized(const EdgeId &edge1, const EdgeId &edge2) const {
            return static_cast <double> (IntersectionSize(edge1, edge2)) / static_cast <double> (g_.length(edge2) + norm_len_);
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

        void ReadEntry (ifstream& fin, const EdgeId& edge, const string& which) override {
            barcode_set_t entry;
            fin >> entry;
            InsertSet(entry, edge, which);
        }

        void WriteEntry (ofstream& fout, const EdgeId& edge, const string& which_end) override {
            VERIFY(which_end == "head" || which_end == "tail");
            barcode_set_t set;
            if (which_end == "head")
                set = GetSetHeads(edge);
            else
                set = GetSetTails(edge);
            fout << edge << ' ' << set << std::endl;
        }

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
        DECL_LOGGER("BarcodeMapperBitset")
    };


    class TrimmableBarcodeMapper : public BarcodeMapper<std::unordered_map <BarcodeId, size_t> > {
        typedef string BarcodeId; //TODO: Encode strings
        typedef std::unordered_map <BarcodeId, size_t> barcode_distribution_t;
        typedef std::unordered_map <EdgeId, barcode_distribution_t> barcode_map_t;
    private:
        using BarcodeMapper::g_;
        using BarcodeMapper::type_;
        barcode_map_t barcode_map_heads;
        barcode_map_t barcode_map_tails;
        size_t tail_threshold_;
        size_t norm_len_;
    public:
        TrimmableBarcodeMapper(const Graph& g, MapperType type,
                               size_t tail_threshold = 10000, size_t norm_len = 10000) :
                BarcodeMapper(g, type), tail_threshold_(tail_threshold), norm_len_(norm_len)
        {
            barcode_map_heads = barcode_map_t();
            barcode_map_tails = barcode_map_t();
        }

        void FillMap(const string &reads_filename, const Index& index,
                     const KmerSubs& kmer_mapper, bool debug_mode = false) override {
            auto lib_vec = GetLibrary(reads_filename);
            auto mapper = std::make_shared<debruijn_graph::NewExtendedSequenceMapper<Graph, Index> >
                    (g_, index, kmer_mapper);

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
                        for (size_t i = 0; i < path.size(); i++) {
                            if (is_at_edge_head(path[i].second)) {
                                InsertBarcode(barcode, path[i].first, "head");
                            }
                            if (is_at_edge_tail(path[i].first, path[i].second))
                                InsertBarcode(barcode, path[i].first, "tail");
                        }
                    }
                    if (debug_mode) {
                        debug_counter++;
                    }
                }
            }
        }

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


        std::pair <double, double> AverageBarcodeCoverage() {
            edge_it_helper helper(g_);
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

        //Delete low abundant barcodes
        void Trim(size_t trimming_threshold) {
            vector <barcode_map_t> maps {barcode_map_heads, barcode_map_tails};
            for (auto map: maps) {
                for (auto entry: map) {
                    barcode_distribution_t& current_distr = map.at(entry.first);
                    for (auto it = current_distr.cbegin(); it != current_distr.cend() ;) {
                        if (it->second < trimming_threshold) {
                            current_distr.erase(it++);
                        }
                        else {
                            ++it;
                        }
                    }
                }
            }
        }

        void SerializeOverallDistribution(std::ofstream& fout) {
            std::unordered_map <size_t, size_t> overall_distr;
            vector <barcode_map_t> maps {barcode_map_heads, barcode_map_tails};
            for (auto map: maps) {
                for (auto entry: map) {
                    barcode_distribution_t& current_distr = map.at(entry.first);
                    for (auto it = current_distr.cbegin(); it != current_distr.cend() ;) {
                        overall_distr[it->second]++;
                    }
                }
            }
            for (auto entry : overall_distr) {
                fout << entry.first << ": " << entry.second << endl;
            }
        }

        void ReadEntry (ifstream& fin, const EdgeId& edge, const string& which) override {
            barcode_distribution_t distr;
            size_t distr_size;
            fin >> distr_size;
            for (size_t i = 0; i < distr_size; ++i) {
                BarcodeId bid;
                size_t abundance;
                fin >> bid >> abundance;
                distr.insert({bid, abundance});
            };
            if (which == "head")
                barcode_map_heads.insert({edge, distr});
            else
                barcode_map_tails.insert({edge, distr});
        }

        void WriteEntry (ofstream& fout, const EdgeId& edge, const string& which) override {
            barcode_distribution_t distr;
            if (which == "head") {
                distr = GetSetHeads(edge);
            }
            else
                distr = GetSetTails(edge);
            size_t distr_size;
            fout << distr.size();
            for (auto entry : distr) {
                fout << entry.first << ' ' << entry.second << endl;
            }
        }

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

        bool is_at_edge_tail(const EdgeId& edge, const omnigraph::MappingRange& range) {
            return range.mapped_range.start_pos + tail_threshold_ > g_.length(edge);
        }

        bool is_at_edge_head(const omnigraph::MappingRange& range) {
            return range.mapped_range.end_pos < tail_threshold_;
        }

        barcode_distribution_t GetSetHeads(const EdgeId &edge) {
            return barcode_map_heads.at(edge);
        }

        barcode_distribution_t GetSetTails(const EdgeId& edge) {
            return barcode_map_tails.at(edge);
        }

        barcode_distribution_t GetSetHeads(const EdgeId &edge) const {
            return barcode_map_heads.at(edge);
        }

        barcode_distribution_t GetSetTails(const EdgeId& edge) const {
            return barcode_map_tails.at(edge);
        }

        void InsertBarcode(const BarcodeId& barcode, const EdgeId& edge, const std::string& which) {
            VERIFY(which == "head" || which == "tail");
            if (which == "head") {
                if (barcode_map_heads.at(edge).find(barcode) == barcode_map_heads[edge].end())
                    barcode_map_heads.at(edge).insert({barcode, 1});
                else
                    barcode_map_heads.at(edge).at(barcode)++;
            }
            else {
                if (barcode_map_tails.at(edge).find(barcode) == barcode_map_tails[edge].end())
                    barcode_map_tails.at(edge).insert({barcode, 1});
                else
                    barcode_map_tails.at(edge).at(barcode)++;
            }
        }

        void InsertSet (barcode_distribution_t& set, const EdgeId& edge, const std::string& which) {
            VERIFY(which == "head" || which == "tail");
            VERIFY(set.size() == max_barcodes);
            if (which == "head")
                barcode_map_heads[edge] = set;
            else
                barcode_map_tails[edge] = set;
        }

        DECL_LOGGER("BarcodeMapperTrimmable")
    };


} //tslr_resolver


