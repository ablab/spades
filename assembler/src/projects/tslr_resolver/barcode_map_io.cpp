#include "barcode_map_io.hpp"

namespace debruijn_graph {
    namespace graphio {
        void SerializeBarcodeMapEntry(size_t e1, const std::unordered_set<BarcodeId>& set, ofstream& file) {
            file << e1 << ' ' << set.size() << std::endl;
            for (auto it : set) {
                file << it << ' ';
            }
            file << std::endl;
        }

        void SerializeMapper(const string& file_name, const tslr_resolver::BarcodeMapper& barcodeMapper) {
            ofstream file;
            file.open(file_name + ".bmap");
            DEBUG("Saving barcode information, " << file_name <<" created");
            VERIFY(file != NULL);
            file << barcodeMapper.size() << std::endl;
            for (auto it = barcodeMapper.cbegin(); it != barcodeMapper.cend(); ++it) {
                SerializeBarcodeMapEntry(it -> first.int_id(), it -> second, file);
            }
        }

        void DeserializeBarcodeMapEntry(ifstream& file, const std::map <size_t, EdgeId>& edge_map, 
                        tslr_resolver::BarcodeMapper& barcodeMapper) {
            size_t edge_id;
            size_t entry_size;

            file >> edge_id;
            file >> entry_size;
            auto edge = edge_map.at(edge_id);
            for (size_t i = 0; i < entry_size; ++i) {
                BarcodeId barcode;
                file >> barcode;
                barcodeMapper.InsertBarcode(barcode, edge);
            }
        }

        void DeserializeMapper(const string& file_name, const std::map <size_t, EdgeId>& edge_map,
                               tslr_resolver::BarcodeMapper& barcodeMapper) {
            ifstream file;
            file.open(file_name);
            DEBUG("Loading barcode information from " << file_name);
            VERIFY(file != NULL);
            size_t map_size;
            file >> map_size;
            for (size_t i = 0; i < map_size; ++i) {
                DeserializeBarcodeMapEntry(file, edge_map, barcodeMapper);
            }
        }



    } //graphio
} //debruijn_graph