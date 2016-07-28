#include <vector>
#include "barcode_mapper.hpp"

namespace debruijn_graph {
    namespace graphio {
        typedef string BarcodeId;
        typedef tslr_resolver::MapperType mtype;
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

        inline string EncodeType (mtype type) {
            switch(type) {
                case mtype::Bitset : return "Bitset";
                case mtype::Trimmable : return "Trimmable";
            }
        }

        inline mtype DecodeType (const string& code) {
            if (code == "Bitset") {
                return mtype::Bitset;
            }
            if (code == "Trimmable") {
                return mtype::Trimmable;
            }
        }

        inline void MapperChooser(std::shared_ptr<BMapper> mapper,
                           mtype map_type, Graph &g) { //TODO Make acceptable reflection
            switch (map_type) {
                case mtype::Bitset :
                    mapper = make_shared<tslr_resolver::BitSetBarcodeMapper> (g, map_type);
                    return;
                case mtype::Trimmable:
                    mapper = make_shared<tslr_resolver::TrimmableBarcodeMapper> (g, map_type);
                    return;
            }
        }

        inline void SerializeMapper(const string& path, shared_ptr<BMapper> barcodeMapper) {
            ofstream file;
            const string file_name = path + ".bmap";
            file.open(file_name);
            DEBUG("Saving barcode information, " << file_name <<" created");
            file << EncodeType(barcodeMapper -> type_);
            file << barcodeMapper->size() << std::endl;
            for (auto it = barcodeMapper->cbegin_heads(); it != barcodeMapper->cend_heads(); ++it) {
                barcodeMapper->WriteEntry(file, it -> first, "head");
            }
            
            file << barcodeMapper->size() << std::endl;
            for (auto it = barcodeMapper->cbegin_tails(); it != barcodeMapper->cend_tails(); ++it) {
                barcodeMapper->WriteEntry(file, it -> first, "tail");
            }
        }

        inline void DeserializeBarcodeMapEntry(ifstream& file, const std::unordered_map <size_t, EdgeId>& edge_map, 
                        shared_ptr<BMapper> barcodeMapper, const std::string& which_end) {
            VERIFY(which_end == "head" || which_end == "tail");
            size_t edge_id;
            file >> edge_id;
            barcodeMapper->ReadEntry(file, edge_map.find(edge_id) -> second, which_end);
        }

        template <class Graph> 
        void DeserializeMapper(const string& path, const std::unordered_map <size_t, EdgeId>& edge_map,
                               shared_ptr<BMapper> barcodeMapper, Graph& g)  {
            ifstream file;
            string file_name = path + ".bmap";
            file.open(file_name);
            INFO("Loading barcode information from " << file_name);
            VERIFY(file != NULL);
            string map_type;
            file >> map_type;
            MapperChooser(barcodeMapper, DecodeType(map_type), g);
            size_t map_size;
            file >> map_size;
            barcodeMapper->InitialFillMap(g);
            for (size_t i = 0; i < map_size; ++i) {
                DeserializeBarcodeMapEntry(file, edge_map, barcodeMapper, "head");
            }
            file >> map_size;
            for (size_t i = 0; i < map_size; ++i) {
                DeserializeBarcodeMapEntry(file, edge_map, barcodeMapper, "tail");
            }
            INFO(barcodeMapper->size());
        }


    } //graphio
} //debruijn_graph
