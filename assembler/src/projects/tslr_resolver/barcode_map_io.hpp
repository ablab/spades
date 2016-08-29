#pragma once

#include <vector>
#include "barcode_mapper.hpp"

namespace debruijn_graph {
    namespace graphio {
        typedef tslr_resolver::BarcodeMapper BMapper;

        template <class Graph>
        std::unordered_map <size_t, typename Graph::EdgeId> MakeEdgeMap(Graph& g) {
            omnigraph::IterationHelper <Graph, typename Graph::EdgeId> helper(g);
            std::unordered_map <size_t, typename Graph::EdgeId> edge_id_map;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                edge_id_map[it -> int_id()] = *it;
            }
            return edge_id_map;
        };

        template <class Graph>
        inline void SerializeMapper(const string& path, const shared_ptr<BMapper>& barcodeMapper, const Graph& g) {
            ofstream file;
            const string file_name = path + ".bmap";
            file.open(file_name);
            DEBUG("Saving barcode information, " << file_name <<" created");
            if (!barcodeMapper) {
                return;
            }
            file << barcodeMapper->size() << std::endl;
            omnigraph::IterationHelper <Graph, typename Graph::EdgeId> helper(g);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                barcodeMapper->WriteEntry(file, *it);
            }
        }

        inline void DeserializeBarcodeMapEntry(ifstream& file, const std::unordered_map <size_t, EdgeId>& edge_map, 
                        shared_ptr<BMapper>& barcodeMapper) {
            size_t edge_id;
            file >> edge_id;
            barcodeMapper->ReadEntry(file, edge_map.find(edge_id) -> second);
        }

        inline void DeserializeMapper(const string& path, const std::unordered_map <size_t, EdgeId>& edge_map,
                               shared_ptr<BMapper>& barcodeMapper)  {
            ifstream file;
            string file_name = path + ".bmap";
            file.open(file_name);
            INFO("Loading barcode information from " << file_name);
            VERIFY(file != NULL);
            size_t map_size;
            file >> map_size;
            for (size_t i = 0; i < map_size; ++i) {
                DeserializeBarcodeMapEntry(file, edge_map, barcodeMapper);
            }
        }


    } //graphio
} //debruijn_graph
