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

        void SerializeMapper(const string& path, const tslr_resolver::BarcodeMapper& barcodeMapper) {
            ofstream file;
            const string file_name = path + ".bmap";
            file.open(file_name);
            DEBUG("Saving barcode information, " << file_name <<" created");
            file << barcodeMapper.size("head") << std::endl;
            for (auto it = barcodeMapper.cbegin_heads(); it != barcodeMapper.cend_heads(); ++it) {
                SerializeBarcodeMapEntry(it -> first.int_id(), it -> second, file);
            }
            
            file << barcodeMapper.size("tail") << std::endl;
            for (auto it = barcodeMapper.cbegin_tails(); it != barcodeMapper.cend_tails(); ++it) {
                SerializeBarcodeMapEntry(it -> first.int_id(), it -> second, file);
            }
        }

        void DeserializeBarcodeMapEntry(ifstream& file, const std::unordered_map <size_t, EdgeId>& edge_map, 
                        tslr_resolver::BarcodeMapper& barcodeMapper, const std::string& which) {
            VERIFY(which == "head" || which == "tail");
            size_t edge_id;
            size_t entry_size;

            file >> edge_id;
            file >> entry_size;
            auto edge = edge_map.at(edge_id);
            for (size_t i = 0; i < entry_size; ++i) {
                BarcodeId barcode;
                file >> barcode;
                barcodeMapper.InsertBarcode(barcode, edge, which);
            }
        }
    } //graphio
} //debruijn_graph
