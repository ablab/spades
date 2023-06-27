#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include <common/barcode_index/barcode_info_extractor.hpp>
#include <iomanip>

namespace online_visualization {

    class BarcodePrinter {
        barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor_;
        const DebruijnEnvironment& curr_env_;
        const string basename_;
        const size_t edge_id_width_;
        const size_t barcode_id_width_;
        const size_t number_of_reads_width_;
        typedef barcode_index::BarcodeId BarcodeId;

    public:
        BarcodePrinter(const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor,
                       const DebruijnEnvironment& curr_env, const string basename) :
                barcode_extractor_(barcode_extractor), curr_env_(curr_env), basename_(basename),
                edge_id_width_(12), barcode_id_width_(9), number_of_reads_width_(5) {}

        void PrintBarcodesOnEdges(const vector<EdgeId>& edges) const {
            stringstream namestream;
            namestream << curr_env_.folder() << "/" << basename_;
            for (const auto& edge: edges) {
                namestream << "_" << edge.int_id();
            }
            INFO("Executing command on " << edges.size() << " edges");
            INFO("Writing barcodes to " << namestream.str());
            ofstream fout(namestream.str());
            vector<BarcodeId> first_last = GetFirstLastIntersection(edges);
            vector<BarcodeId> barcode_intersection = GetBarcodeIntersection(edges);
            vector<BarcodeId> barcode_union = GetBarcodeUnion(edges);
            fout << "Edges: ";
            for (const auto& edge: edges) {
                fout << edge.int_id() << ", ";
            }
            fout << endl;
            fout << "Intersection of first with last, size: " << first_last.size() << endl;
            PrintRows(first_last, edges, fout);
            fout << "Intersection, size: " << barcode_intersection.size() << endl;
            PrintRows(barcode_intersection, edges, fout);
            fout << "Union, size: " << barcode_union.size() << endl;
            PrintRows(barcode_union, edges, fout);
        }



    private:
        vector <BarcodeId> GetBarcodeIntersection(const vector<EdgeId>& edges) const {
            VERIFY(edges.size() > 1);
            vector<BarcodeId> current_intersection = barcode_extractor_.GetBarcodes(edges[0]);
            for (auto it = std::next(edges.begin()); it != edges.end(); ++it) {
                vector <BarcodeId> new_intersection;
                vector <BarcodeId> edge_barcodes = barcode_extractor_.GetBarcodes(*it);
                std::set_intersection(current_intersection.begin(), current_intersection.end(),
                                      edge_barcodes.begin(), edge_barcodes.end(),
                                      std::back_inserter(new_intersection));
                current_intersection.clear();
                std::move(new_intersection.begin(), new_intersection.end(), std::back_inserter(current_intersection));
            }
            return current_intersection;
        }

        vector <BarcodeId> GetBarcodeUnion(const vector<EdgeId>& edges) const {
            VERIFY(edges.size() > 1);
            vector<BarcodeId> current_union = barcode_extractor_.GetBarcodes(edges[0]);
            for (auto it = std::next(edges.begin()); it != edges.end(); ++it) {
                vector <BarcodeId> new_union;
                vector <BarcodeId> edge_barcodes = barcode_extractor_.GetBarcodes(*it);
                std::set_union(current_union.begin(), current_union.end(),
                               edge_barcodes.begin(), edge_barcodes.end(),
                               std::back_inserter(new_union));
                current_union.clear();
                std::move(new_union.begin(), new_union.end(), std::back_inserter(current_union));
            }
            return current_union;
        }

        vector <BarcodeId> GetFirstLastIntersection(const vector<EdgeId>& edges) const {
            VERIFY(edges.size() > 1);
            EdgeId first = edges[0];
            EdgeId last = edges.back();
            return barcode_extractor_.GetSharedBarcodes(first, last);
        }

        void PrintRows(const vector<BarcodeId> barcodes, const vector<EdgeId>& edges, ofstream& fout) const {
            for (const auto& barcode: barcodes) {
                PrintRow(barcode, edges, fout);
            }
        }

        void PrintRow(const BarcodeId& barcode, const vector<EdgeId>& edges, ofstream& fout) const {
            PrintElement(barcode, barcode_id_width_, fout);
            for (const auto& edge: edges) {
                if (not barcode_extractor_.HasBarcode(edge, barcode)) {
                    const size_t reads = 0;
                    const auto bitset_size = barcode_extractor_.GetNumberOfBins(edge);
                    boost::dynamic_bitset<> bitset;
                    bitset.resize(bitset_size);
                    const size_t bitset_width = bitset_size + 5;
                    PrintElement(reads, number_of_reads_width_, fout);
                    PrintElement(bitset, bitset_width, fout);
                } else {
                    const size_t reads = barcode_extractor_.GetNumberOfReads(edge, barcode);
                    PrintElement(reads, number_of_reads_width_, fout);
                }
            }
            fout << endl;
        }

        template <class T> void PrintElement(const T& element, const size_t width, ofstream& fout) const  {
            fout << std::left << std::setw(static_cast<int>(width)) << element;
        }
    };

    class PrintBarcodesCommand : public LocalCommand<DebruijnEnvironment> {
    private:


    protected:
        size_t MinArgNumber() const {
            return 2;
        }

        bool CheckCorrectness(const vector<string>& args) const {
            return CheckEnoughArguments(args);
        }

    public:
        string Usage() const {
            string answer;
            answer = answer + "Command `paths` \n" +
                     "Usage:\n" +
                     "> print_barcodes <edgeid1> <edgeid2> ... <edgeidn> \n" +
                     " This command prints barcodes aligned to every edge from the list along with compressed read cloud alignments.";
            return answer;
        }

        PrintBarcodesCommand() : LocalCommand<DebruijnEnvironment>("print_barcodes")
        {
        }

        void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
            const vector<string>& args = arg_list.GetAllArguments();
            if (!CheckCorrectness(args))
                return;
            TRACE("Executing `print_barcodes` command");
            vector<EdgeId> edges;
            for (size_t i = 1; i < args.size(); ++i) {
                size_t edge_id = GetInt(args[i]);
                if (!CheckEdgeExists(curr_env.finder(), edge_id)) {
                    INFO("Edge " << edge_id << " was not found in graph.");
                    return;
                }
                edges.push_back(curr_env.finder().ReturnEdgeId(edge_id));
            }
            string basename = "barcode_table_";
            make_dir(curr_env.folder());
            barcode_index::FrameBarcodeIndexInfoExtractor extractor(curr_env.GetBarcodeExtractor(), curr_env.graph());
            BarcodePrinter printer(extractor, curr_env, basename);
            printer.PrintBarcodesOnEdges(edges);
        }

        DECL_LOGGER("PrintBarcodes");
    };

}
