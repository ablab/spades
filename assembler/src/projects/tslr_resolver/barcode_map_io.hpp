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

        void SerializeBarcodeMapEntry(size_t e1, const std::unordered_set<BarcodeId>& set, ofstream& file);

        void SerializeMapper(const string& path, const tslr_resolver::BarcodeMapper& barcodeMapper);

        void DeserializeBarcodeMapEntry(ifstream& file, const std::unordered_map <size_t, EdgeId>& edge_map, 
                        tslr_resolver::BarcodeMapper& barcodeMapper, const std::string& which);

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
            INFO(barcodeMapper.size("head"));
            INFO(barcodeMapper.size("tail"));
        }


    } //graphio
} //debruijn_graph
