#include <vector>
#include "barcode_mapper.hpp"

namespace debruijn_graph {
    namespace graphio {
        typedef string BarcodeId;

        template <class Graph>
        std::unordered_map <size_t, typename Graph::EdgeId> MakeEdgeMap(Graph& g) {
            omnigraph::IterationHelper <Graph, typename Graph::EdgeId> helper(g);
            std::unordered_map <size_t, typename Graph::EdgeId> edge_id_map;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                edge_id_map[it -> int_id()] = *it;
            }
            return edge_id_map;
        };

        inline void SerializeBarcodeMapEntry(size_t e1, const tslr_resolver::barcode_set_t& set, ofstream& file) {
            file << e1 << ' ' << set << std::endl;
            file << std::endl;
        }

        inline void SerializeMapper(const string& path, const tslr_resolver::BarcodeMapper& barcodeMapper) {
            ofstream file;
            const string file_name = path + ".bmap";
            file.open(file_name);
            DEBUG("Saving barcode information, " << file_name <<" created");
            file << barcodeMapper.size() << std::endl;
            for (auto it = barcodeMapper.cbegin_heads(); it != barcodeMapper.cend_heads(); ++it) {
                SerializeBarcodeMapEntry(it -> first.int_id(), it -> second, file);
            }
            
            file << barcodeMapper.size() << std::endl;
            for (auto it = barcodeMapper.cbegin_tails(); it != barcodeMapper.cend_tails(); ++it) {
                SerializeBarcodeMapEntry(it -> first.int_id(), it -> second, file);
            }
        }

        inline void DeserializeBarcodeMapEntry(ifstream& file, const std::unordered_map <size_t, EdgeId>& edge_map, 
                        tslr_resolver::BarcodeMapper& barcodeMapper, const std::string& which) {
            VERIFY(which == "head" || which == "tail");
            size_t edge_id;
            tslr_resolver::barcode_set_t entry;

            file >> edge_id;
            file >> entry;
            auto edge = edge_map.at(edge_id);
            barcodeMapper.InsertSet(entry, edge, which);
        }

        template <class Graph> 
        void DeserializeMapper(const string& path, const std::unordered_map <size_t, EdgeId>& edge_map,
                               tslr_resolver::BarcodeMapper& barcodeMapper, Graph& g)  {
            ifstream file;
            string file_name = path + ".bmap";  //TODO: Get Stage name somehow
            file.open(file_name);
            INFO("Loading barcode information from " << file_name);
            VERIFY(file != NULL);
            size_t map_size;
            file >> map_size;
            barcodeMapper.InitialFillMap(g);
            for (size_t i = 0; i < map_size; ++i) {
                DeserializeBarcodeMapEntry(file, edge_map, barcodeMapper, "head");
            }
            file >> map_size;
            for (size_t i = 0; i < map_size; ++i) {
                DeserializeBarcodeMapEntry(file, edge_map, barcodeMapper, "tail");
            }
            INFO(barcodeMapper.size());
            INFO(barcodeMapper.size());
        }


    } //graphio
} //debruijn_graph
