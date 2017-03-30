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

        void PrintBarcodesOnEdges(const EdgeId& first_edge, const EdgeId& second_edge) const {
            stringstream namestream;
            namestream << curr_env_.folder() << "/" << basename_ << "_" << first_edge.int_id() << "_" << second_edge.int_id();
            INFO("Writing barcodes to " << namestream.str());
            ofstream fout(namestream.str());
            vector <BarcodeId> barcode_intersection = barcode_extractor_.GetSharedBarcodes(first_edge, second_edge);
            vector <BarcodeId> first_barcodes = barcode_extractor_.GetBarcodes(first_edge);
            vector <BarcodeId> second_barcodes = barcode_extractor_.GetBarcodes(second_edge);
            vector <BarcodeId> barcode_union;
            std::set_union(first_barcodes.begin(), first_barcodes.end(),
                           second_barcodes.begin(), second_barcodes.end(),
                           std::back_inserter(barcode_union));
            fout << "Edges " << first_edge.int_id() << " and " << second_edge.int_id() << endl;
            fout << "Intersection, size: " << barcode_intersection.size() << endl;
            PrintRows(barcode_intersection, first_edge, second_edge, fout);
//            fout << "Union, size: " << barcode_union.size() << endl;
//            PrintRows(barcode_union, first_edge, second_edge, fout);
        }

    private:
        void PrintRows(const vector<BarcodeId> barcodes, const EdgeId& first_edge,
                       const EdgeId& second_edge, ofstream& fout) const {
            for (const auto& barcode: barcodes) {
                PrintRow(barcode, first_edge, second_edge, fout);
            }
        }

        void PrintRow(const BarcodeId& barcode, const EdgeId& first_edge, const EdgeId& second_edge, ofstream& fout) const {
            const size_t first_reads = barcode_extractor_.GetNumberOfReads(first_edge, barcode);
            const size_t second_reads = barcode_extractor_.GetNumberOfReads(second_edge, barcode);
            const auto first_bitset = barcode_extractor_.GetBitSet(first_edge, barcode);
            const auto second_bitset = barcode_extractor_.GetBitSet(second_edge, barcode);
            const size_t first_bitset_width = first_bitset.size() + 5;
            const size_t second_bitset_width = second_bitset.size() + 5;
            PrintElement(barcode, barcode_id_width_, fout);
            PrintElement(first_reads, number_of_reads_width_, fout);
            PrintElement(first_bitset, first_bitset_width, fout);
            PrintElement(second_reads, number_of_reads_width_, fout);
            PrintElement(second_bitset, second_bitset_width, fout);
            fout << endl;
        }

        template <class T> void PrintElement(const T& element, const size_t width, ofstream& fout) const  {
            fout << std::left << std::setw(static_cast<int>(width)) << element;
        }
    };

    class PrintBarcodes : public LocalCommand<DebruijnEnvironment> {
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
                     "> print_barcodes <edgeid1> <edgeid2> \n" +
                     " This command prints barcodes aligned to edge1 and edge2 along with compressed alignments.";
            return answer;
        }

        PrintBarcodes() : LocalCommand<DebruijnEnvironment>("print_barcodes")
        {
        }

        void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
            const vector<string>& args = arg_list.GetAllArguments();
            if (!CheckCorrectness(args))
                return;
            TRACE("Executing `print_barcodes` command");
            size_t edgeid1 = GetInt(args[1]);
            size_t edgeid2 = GetInt(args[2]);
            if (!CheckEdgeExists(curr_env.finder(), edgeid1) or !CheckEdgeExists(curr_env.finder(), edgeid2))
                return;
            EdgeId edge1 = curr_env.finder().ReturnEdgeId(edgeid1);
            EdgeId edge2 = curr_env.finder().ReturnEdgeId(edgeid2);
            string basename = "barcode_table_";
            make_dir(curr_env.folder());
            barcode_index::FrameBarcodeIndexInfoExtractor extractor(curr_env.GetBarcodeExtractor(), curr_env.graph());
            BarcodePrinter printer(extractor, curr_env, basename);
            printer.PrintBarcodesOnEdges(edge1, edge2);
        }

        DECL_LOGGER("PrintBarcodes");
    };

}
