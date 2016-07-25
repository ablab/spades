#include "barcode_map_io.hpp"

namespace debruijn_graph {
    namespace graphio {
        void SerializeBarcodeMapEntry(size_t e1, const tslr_resolver::barcode_set_t& set, ofstream& file) {
            file << e1 << ' ' << set << std::endl;
            file << std::endl;
        }

        void SerializeMapper(const string& path, const tslr_resolver::BarcodeMapper& barcodeMapper) {
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

        void DeserializeBarcodeMapEntry(ifstream& file, const std::unordered_map <size_t, EdgeId>& edge_map, 
                        tslr_resolver::BarcodeMapper& barcodeMapper, const std::string& which) {
            VERIFY(which == "head" || which == "tail");
            size_t edge_id;
            tslr_resolver::barcode_set_t entry;

            file >> edge_id;
            file >> entry;
            auto edge = edge_map.at(edge_id);
            barcodeMapper.InsertSet(entry, edge, which);
        }
    } //graphio
} //debruijn_graph
