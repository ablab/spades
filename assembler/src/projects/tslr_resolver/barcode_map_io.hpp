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

        void SerializeMapper(const string& file_name, const tslr_resolver::BarcodeMapper& barcodeMapper);

        void DeserializeBarcodeMapEntry(ifstream& file, const std::unordered_map <size_t, EdgeId>& edge_map, 
                        tslr_resolver::BarcodeMapper& barcodeMapper);

        void DeserializeMapper(const string& file_name, const std::unordered_map <size_t, EdgeId>& edge_map,
                               tslr_resolver::BarcodeMapper& barcodeMapper);

    } //graphio
} //debruijn_graph
