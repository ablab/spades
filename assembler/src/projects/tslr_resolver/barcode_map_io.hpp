#pragma once

#include <vector>
#include "barcode_mapper.hpp"

using namespace tslr_resolver;

namespace debruijn_graph {
    namespace graphio {

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
        inline void SerializeMapper(const string& path, const shared_ptr<BarcodeMapper>& barcodeMapper, const Graph& g) {
            ofstream file;
            const string file_name = path + ".bmap";
            file.open(file_name);
            if (!barcodeMapper || barcodeMapper->IsEmpty()) {
                return;
            }
            INFO(barcodeMapper->size())
            file << barcodeMapper->size() << std::endl;
            omnigraph::IterationHelper <Graph, typename Graph::EdgeId> helper(g);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                barcodeMapper->WriteEntry(file, *it);
            }
        }

        inline void DeserializeBarcodeMapEntry(ifstream& file, const std::unordered_map <size_t, EdgeId>& edge_map, 
                        shared_ptr<BarcodeMapper>& barcodeMapper) {
            size_t edge_id;
            file >> edge_id;
            barcodeMapper->ReadEntry(file, edge_map.find(edge_id) -> second);
        }

        template <class Graph>
        void DeserializeMapper(const string& path, const std::unordered_map <size_t, EdgeId>& edge_map,
                               shared_ptr<BarcodeMapper>& barcodeMapper, Graph& g)  {
            ifstream file;
            string file_name = path + ".bmap";
            file.open(file_name);
            INFO("Loading barcode information from " << file_name)
            VERIFY(file != NULL);
            HeadTailMapperBuilder<SimpleBarcodeEntry> mapper_builder(g, cfg::get().ts_res.edge_tail_len);
            INFO("Built mapper")
            barcodeMapper = mapper_builder.GetMapper();
            INFO(barcodeMapper->size())
            if(file.peek() == std::ifstream::traits_type::eof()) {
                return;
            }
            size_t map_size;
            file >> map_size;
            for (size_t i = 0; i < map_size; ++i) {
                DeserializeBarcodeMapEntry(file, edge_map, barcodeMapper);
            }
        }


    } //graphio
} //debruijn_graph
