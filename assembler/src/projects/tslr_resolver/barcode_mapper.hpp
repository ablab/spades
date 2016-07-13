#pragma once

#include <boost/unordered>
#include <memory>
#include <utility>
#include <fstream>
#include <string>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include "io/reads/paired_read.hpp"

using std::string;

namespace tslr_resolver {
    typedef debruijn_graph::conj_graph_pack graph_pack;
    typedef debruijn_graph::ConjugateDeBruijnGraph Graph;
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;
    typedef string BarcodeId;
    typedef boost::unordered_set <BarcodeId> BarcodeSet;
    typedef boost::unordered_map <EdgeId, BarcodeSet> barcode_map_t;
    typedef omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper;


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

    template <typename ReadType>
    class barcode_mapper {
    private:
        barcode_map_t barcode_map_;
        string reads_filename_;
        const graph_pack &gp_;
    public:
        barcode_mapper(const graph_pack &gp, const string &reads_filename)
        try : reads_filename_(reads_filename), gp_(gp) {
            ConstructMap();
        }
        catch (std::exception const& e) {
            std::cerr << "Exception caught " << e.what() << std::endl;
        }


        BarcodeSet GetSet(const EdgeId &edge) const {
            return barcode_map_.at(edge);
        }

        size_t IntersectionSize(const EdgeId &edge1, const EdgeId &edge2) const {
            size_t ans = 0;
            auto Set1 = GetSet(edge1);
            auto Set2 = GetSet(edge2);
            for (auto it = Set1.begin(); it != Set1.end(); ++it) {
                auto it2 = Set2.find(*it);
                if (it2 != Set2.end()) {
                    ans++;
                }
            }
            return ans;
        }

    private:
        void ConstructMap() {
            std::ifstream fin;
            fin.open(reads_filename_);
            std::vector <barcode_library> lib_vec;
            while (!fin.eof()) {
                barcode_library lib;
                fin >> lib.barcode_;
                fin >> lib.left_;
                fin >> lib.right_;
                lib_vec.push_back(lib);
            }

            edge_it_helper helper(gp_.g);
            barcode_map_ = barcode_map_t();
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                BarcodeSet set;
                barcode_map_.insert({*it, set});
            }
            auto mapper = debruijn_graph::MapperInstance <graph_pack> (gp_);

            for (auto lib: lib_vec) {
                std::string barcode = lib.barcode_;
                auto paired_read_stream = io::SeparatePairedReadStream(lib.left_, lib.right_, 0);
                io::PairedRead read;
                while (!paired_read_stream.eof()) {
                    paired_read_stream >> read;
                    auto path_first = mapper -> MapRead(read.first());
                    auto path_second = mapper -> MapRead(read.second());
                    for(size_t i = 0; i < path_first.size(); i++) {
                        barcode_map_.at(path_first[i].first).insert(barcode);
                    }
                    for(size_t i = 0; i < path_second.size(); i++) {
                        barcode_map_.at(path_second[i].first).insert(barcode);
                    }
                }
            }
        }
    };


} //tslr_resolver
