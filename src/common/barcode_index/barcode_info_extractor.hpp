//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2017-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_index.hpp"

namespace barcode_index {

    /**
     * BarcodeIndexInfoExtractor extracts information from BarcodeIndex
     */
    template<class Graph, class BarcodeEntryT>
    class BarcodeIndexInfoExtractor {
    public:
        typedef typename BarcodeEntryT::barcode_distribution_t distribution_t;
        typedef typename BarcodeEntryT::barcode_info_t barcode_info_t;
        typedef typename distribution_t::key_type barcode_info_key_t;
        typedef typename distribution_t::mapped_type barcode_info_value_t;
        typedef typename distribution_t::value_type barcode_info_pair_t;
        typedef typename barcode_index::BarcodeIndex<Graph, BarcodeEntryT> BarcodeIndexT;
        typedef typename Graph::EdgeId EdgeId;
        typedef typename omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper;
    protected:
        const BarcodeIndexT &index_;
        const Graph &g_;
    public:
        BarcodeIndexInfoExtractor(const BarcodeIndexT &index, const Graph &g) :
                index_(index), g_(g) {}

        /**
         *
         * @param edge
         * @return Number of barcodes contained by the edge
         */
        size_t GetNumberOfBarcodes(EdgeId edge) const {
            return index_.GetEntry(edge).Size();
        }

        /**
         *
         * @param edge1
         * @param edge2
         * @return List of barcodes shared by edge1 and edge2
         */
        std::vector<BarcodeId> GetSharedBarcodes (EdgeId edge1, EdgeId edge2) const {
            std::vector<BarcodeId> intersection;
            for (auto it = intersection_iterator_begin(edge1, edge2); it != intersection_iterator_end(edge1, edge2); ++it) {
                intersection.push_back((*it).key_);
            }
            return intersection;
        }

        /**
         *
         * @param edge1
         * @param edge2
         * @return Number of barcodes shared by edge1 and edge2
         */
        size_t GetNumberOfSharedBarcodes (EdgeId edge1, EdgeId edge2) const {
            size_t result = 0;
            for (auto it = intersection_iterator_begin(edge1, edge2); it != intersection_iterator_end(edge1, edge2); ++it) {
                ++result;
            }
            return result;
        }

        /**
         * @param edge
         * @param barcode
         * @return True if the edge contains the barcode
         */
        bool HasBarcode(EdgeId edge, BarcodeId barcode) const {
            return index_.GetEntry(edge).has_barcode(barcode);
        }

        /**
         *
         * @return Average number of barcodes contained on the long edges of the graph
         */
        double AverageBarcodeCoverage(size_t length_threshold) const {
            edge_it_helper helper(g_);
            size_t barcodes_overall = 0;
            size_t long_edges = 0;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                if (g_.length(*it) > length_threshold) {
                    long_edges++;
                    barcodes_overall += GetNumberOfBarcodes(*it);
                }
            }
            DEBUG("tails: " + std::to_string(barcodes_overall));
            DEBUG("Long edges" + long_edges);
            return static_cast <double> (barcodes_overall) / static_cast <double> (long_edges);
        }

        double GetIntersectionSizeNormalizedByUnion(EdgeId edge1, EdgeId edge2) const {
            if (GetUnionSize(edge1, edge2)) {
                return static_cast <double> (GetNumberOfSharedBarcodes(edge1, edge2)) /
                       static_cast <double> (GetUnionSize(edge1, edge2));
            }
            return 0;
        }

        virtual double GetIntersectionSizeNormalizedBySecond(EdgeId edge1, EdgeId edge2) const {
            if (GetNumberOfBarcodes(edge2) > 0) {
                return static_cast <double> (GetNumberOfSharedBarcodes(edge1, edge2)) /
                       static_cast <double> (GetNumberOfBarcodes(edge2));
            }
            return 0;
        }

        size_t GetUnionSize(EdgeId edge1, EdgeId edge2) const {
            auto it_tail = index_.GetEntryTailsIterator(edge1);
            auto it_head = index_.GetEntryHeadsIterator(edge2);
            return (it_tail->second).GetUnionSize(it_head->second);
        }

        std::vector<BarcodeId> GetBarcodes(EdgeId edge) const {
            std::vector <BarcodeId> result;
            auto copy_barcode_id = [&result](const barcode_info_pair_t& entry)
                        {result.push_back(entry.first); };
            std::for_each(barcode_iterator_begin(edge), barcode_iterator_end(edge), copy_barcode_id);
            return result;
        }

        typename distribution_t::const_iterator barcode_iterator_begin(EdgeId edge) const {
            auto entry_it = index_.GetEntryHeadsIterator(edge);
            return entry_it->second.begin();
        }

        typename distribution_t::const_iterator barcode_iterator_end(EdgeId edge) const {
            auto entry_it = index_.GetEntryHeadsIterator(edge);
            return entry_it->second.end();
        }

        /**
         * Proxy class representing a pair of references to BarcodeInfo
         */
        struct IntersectionData {
            /**
             * BarcodeId of the shared barcode
             */
            const barcode_info_key_t key_;
            /**
             * Info corresponding to the shared barcode and the first edge
             */
            const barcode_info_value_t& info_first_;

            /**
             * Info corresponding to the shared barcode and the second edge
             */
            const barcode_info_value_t& info_second_;

            IntersectionData(const barcode_info_key_t key,
                             const barcode_info_value_t& info_first,
                             const barcode_info_value_t& info_second) :
                    key_(key), info_first_(info_first), info_second_(info_second) {}
        };

        /**
         * Iterator over shared barcodes of two edges.
         * Dereferencing returns proxy object of type IntersectionData
         * @note Since it is not an iterator over a container there is no -> operator.
         */
        class const_intersection_iterator {
        public:
            typedef IntersectionData value_type;
            typedef IntersectionData reference;
            typedef IntersectionData* pointer;
            typedef int difference_type;
            typedef std::input_iterator_tag iterator_category;
            typedef typename distribution_t::const_iterator entry_iterator;

        public:
            const_intersection_iterator(entry_iterator first, entry_iterator second,
                                        entry_iterator first_end, entry_iterator second_end) :
                    first_(first), second_(second), first_end_(first_end), second_end_(second_end) {}

            //todo optimize with lower bounds?
            const_intersection_iterator operator++() {
                if (first_ == first_end_ and second_ == second_end_) {
                    ++first_;
                    ++second_;
                    return *this;
                }
                if (get_first_key() == get_second_key()) {
                    if (second_ != second_end_) {
                        ++second_;
                    } else {
                        ++first_;
                    }
                }
                while (get_first_key() != get_second_key() and (first_ != first_end_ or second_ != second_end_)) {
                    while (get_first_key() < get_second_key() and first_ != first_end_) {
                        ++first_;
                    }
                    while (get_second_key() < get_first_key() and second_ != second_end_) {
                        ++second_;
                    }
                    if ((first_ == first_end_ and get_second_key() > get_first_key()) or
                            (second_ == second_end_ and get_first_key() > get_second_key())) {
                        first_ = first_end_;
                        second_ = second_end_;
                    }
                }
                if (get_first_key() == get_second_key()) {
                    return *this;
                }
                VERIFY(first_ == first_end_ and second_ == second_end_);
                if (get_first_key() != get_second_key()) {
                    ++first_;
                    ++second_;
                }
                return *this;
            }

            const_intersection_iterator operator++(int) {
                const_intersection_iterator result = *this;
                ++(*this);
                return result;
            }

            reference operator* () {
                VERIFY(get_first_key() == get_second_key());
                return IntersectionData(get_first_key(), first_->second, second_->second);
            }

            bool operator== (const const_intersection_iterator& other) {
                return first_ == other.first_ and second_ == other.second_ and
                                                  first_end_ == other.first_end_ and second_end_ == other.second_end_;
            }

            bool operator!= (const const_intersection_iterator& other) {
                return *this != other;
            }

        private:
            entry_iterator first_;
            entry_iterator second_;
            entry_iterator first_end_;
            entry_iterator second_end_;

            inline barcode_info_key_t get_first_key() {
                return first_->first;
            }
            inline barcode_info_key_t get_second_key() {
                return second_->first;
            }
        };

        //todo remove end decrement?
        const_intersection_iterator intersection_iterator_begin(EdgeId first, EdgeId second) const {
            if (GetNumberOfBarcodes(first) == 0 or GetNumberOfBarcodes(second) == 0) {
                return intersection_iterator_end(first, second);
            }
            auto first_begin = barcode_iterator_begin(first);
            auto first_end = barcode_iterator_end(first);
            auto second_begin = barcode_iterator_begin(second);
            auto second_end = barcode_iterator_end(second);
            const_intersection_iterator prebegin(first_begin, second_begin, --first_end, --second_end);
            if (barcode_iterator_begin(first)->first == barcode_iterator_begin(second)->first)
                return prebegin;
            return ++prebegin;
        }

        const_intersection_iterator intersection_iterator_end(EdgeId first, EdgeId second) const {
            auto first_end = barcode_iterator_end(first);
            auto second_end = barcode_iterator_end(second);
            if (GetNumberOfBarcodes(first) == 0 or GetNumberOfBarcodes(second) == 0) {
                return const_intersection_iterator(first_end, second_end, first_end, second_end);
            }
            auto first_end_copy = first_end;
            auto second_end_copy = second_end;
            const_intersection_iterator result(first_end_copy, second_end_copy, --first_end, --second_end);
            return result;
        }

    protected:

        /**
         *
         * @param edge
         * @param barcode
         * @return barcode info corresponding to given edge
         */
        const barcode_info_t& GetInfo(EdgeId edge, BarcodeId barcode) const {
            VERIFY(HasBarcode(edge, barcode));
            const BarcodeEntryT& entry = GetEntry(edge);
            return entry.get_barcode(barcode)->second;
        }

        const BarcodeEntryT& GetEntry(EdgeId edge) const {
            return index_.GetEntry(edge);
        }

    };

/**
 * Specialization of BarcodeIndexInfoExtractor for FrameBarcodeInfo type.
 * @see FrameBarcodeInfo
 */
template<class Graph>
class FrameBarcodeIndexInfoExtractorTemplate : public BarcodeIndexInfoExtractor<Graph, FrameEdgeEntry<Graph>> {
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename barcode_index::FrameBarcodeIndex<Graph> FrameBarcodeIndexT;

    FrameBarcodeIndexInfoExtractorTemplate(const FrameBarcodeIndexT &index, const Graph &g) :
            BarcodeIndexInfoExtractor<Graph, FrameEdgeEntry<Graph>>(index, g) {}

    /**
     *
     * @param edge
     * @param barcode
     * @return number of barcoded reads aligned to the edge
     */
    size_t GetNumberOfReads(EdgeId edge, const BarcodeId &barcode) const {
        return this->GetInfo(edge, barcode).GetCount();
    }

    /**
     *
     * @param edge
     * @param barcode
     * @return leftmost barcoded bin of the edge
     */
    size_t GetLeftBin(EdgeId edge, const BarcodeId &barcode) const {
        return this->GetInfo(edge, barcode).GetLeftMost();
    }

    /**
     *
     * @param edge
     * @param barcode
     * @return rightmost barcoded bin of the edge
     */
    size_t GetRightBin(EdgeId edge, const BarcodeId &barcode) const {
        return this->GetInfo(edge, barcode).GetRightMost();
    }

    /**
     * @param edge
     * @return length of the bin
     */
    size_t GetBinLength(EdgeId edge) const {
        return this->GetEntry(edge).GetFrameSize();
    }

    /**
     *
     * @param edge
     * @return number of bins on the edge
     */
    size_t GetNumberOfBins(EdgeId edge) const {
        return this->GetEntry(edge).GetNumberOfFrames();
    }

    /**
     * @param first first edge
     * @param second second edge
     * @param shared_threshold minimal number of barcodes shared by first and second
     * @param count_threshold edge contains barcode iff there are at least count_threshold reads aligned to the edge
     * @param gap_threshold clouds located at the beginning of the first or at the end of the second edge are discarded.
     *      Cloud is located in the beginning of the edge if it is not aligned to the last gap_threshold nucleotides of the edge.
     * @return true if there are at least shared_threshold barcodes which pass requirements determined by count_threshold and gap_threshold.
     */
    bool AreEnoughSharedBarcodesWithFilter (EdgeId first,
                                            EdgeId second,
                                            size_t shared_threshold,
                                            size_t count_threshold,
                                            size_t gap_threshold) const {
        size_t current = 0;
        for (auto it = intersection_iterator_begin(first, second); it != intersection_iterator_end(first, second); ++it) {
            BarcodeId barcode = (*it).key_;
            bool is_in_the_end_of_first = g_.length(first) <= gap_threshold or
                                                 GetMaxPos(first, barcode) + gap_threshold > g_.length(first);
            bool is_in_the_beginning_of_second = g_.length(second) <= gap_threshold or
                                                 GetMinPos(second, barcode) < gap_threshold;
            bool enough_count = (*it).info_first_.GetCount() >= count_threshold and
                                (*it).info_second_.GetCount() >= count_threshold;
            if (is_in_the_end_of_first and is_in_the_beginning_of_second and enough_count) {
                ++current;
            }
            if (current > shared_threshold) {
                return true;
            }
        }
        return false;
    }

    /**
     * @param first first edge
     * @param second second edge
     * @param count_threshold edge contains barcode iff there are at least count_threshold reads aligned to the edge
     * @param gap_threshold clouds located at the beginning of the first or at the end of the second edge are discarded.
     *      Cloud is located in the beginning of the edge if it is not aligned to the last gap_threshold nucleotides of the edge.
     * @return number of barcodes which pass requirements determined by count_threshold and gap_threshold.
     */
    size_t CountSharedBarcodesWithFilter (EdgeId first,
                                         EdgeId second,
                                         size_t count_threshold,
                                         size_t gap_threshold) const {
        return GetSharedBarcodesWithFilter(first, second, count_threshold, gap_threshold).size();
    }

    std::vector<BarcodeId> GetSharedBarcodesWithFilter(EdgeId first, EdgeId second,
                                                  size_t count_threshold, size_t gap_threshold) const {
        std::vector<BarcodeId> intersection;
        for (auto it = this->intersection_iterator_begin(first, second); it != this->intersection_iterator_end(first, second); ++it) {
            auto barcode = (*it).key_;
            bool is_in_the_end_of_first = g_.length(first) <= gap_threshold or
                GetMaxPos(first, barcode) + gap_threshold > g_.length(first);
            bool is_in_the_beginning_of_second = g_.length(second) <= gap_threshold or
                GetMinPos(second, barcode) < gap_threshold;
            bool enough_count = (*it).info_first_.GetCount() >= count_threshold and
                (*it).info_second_.GetCount() >= count_threshold;
            if (is_in_the_end_of_first and is_in_the_beginning_of_second and enough_count) {
                intersection.push_back(barcode);
            }
        }
        return intersection;
    }

    std::vector<BarcodeId> GetBarcodesFromHead(EdgeId edge, size_t count_threshold, size_t right) const {
        std::vector<BarcodeId> barcodes;
        size_t bin_length = GetBinLength(edge);
        for (auto it = this->barcode_iterator_begin(edge); it != this->barcode_iterator_end(edge); ++it) {
            BarcodeId barcode = it->first;
            size_t left_pos = it->second.GetLeftMost() *  bin_length;
            size_t reads = it->second.GetCount();
            if (left_pos <= right and reads >= count_threshold) {
                barcodes.push_back(barcode);
            }
        }
        return barcodes;
    }

    std::vector<std::pair<BarcodeId, size_t>> GetBarcodesAndCountsFromHead(EdgeId edge,
                                                                           size_t count_threshold,
                                                                           size_t right) const {
        std::vector<std::pair<BarcodeId, size_t>> barcodes;
        size_t bin_length = GetBinLength(edge);
        for (auto it = this->barcode_iterator_begin(edge); it != this->barcode_iterator_end(edge); ++it) {
            BarcodeId barcode = it->first;
            size_t left_pos = it->second.GetLeftMost() *  bin_length;
            size_t reads = it->second.GetCount();
            if (left_pos <= right and reads >= count_threshold) {
                barcodes.emplace_back(barcode, reads);
            }
        }
        return barcodes;
    };

    std::vector<BarcodeId> GetBarcodesFromRange(EdgeId edge, size_t count_threshold,
                                                size_t left, size_t right) const {
        std::vector<BarcodeId> barcodes;
        size_t bin_length = GetBinLength(edge);
        for (auto it = this->barcode_iterator_begin(edge); it != this->barcode_iterator_end(edge); ++it) {
            BarcodeId barcode = it->first;
            size_t left_pos = it->second.GetLeftMost() *  bin_length;
            size_t right_pos = it->second.GetRightMost() * bin_length;
            TRACE("Bin length: " << bin_length);
            TRACE("Left raw: " << left_pos);
            TRACE("Leftmost position: " << left_pos);
            TRACE("Rightmost position: " << right_pos);
            size_t reads = it->second.GetCount();
            TRACE("Reads: " << reads);
            if (left_pos <= right and right_pos >= left and reads >= count_threshold) {
                barcodes.push_back(barcode);
            }
        }
        return barcodes;
    }

    /**
     *
     * @param edge
     * @param barcode
     * @return Estimated first position of the cloud defined by the barcode and the edge (not the first bin, but the first nucleotide)
     */
    size_t GetMinPos(EdgeId edge, BarcodeId barcode) const {
        const FrameEdgeEntry<Graph> &entry = this->GetEntry(edge);
        const FrameBarcodeInfo& info = this->GetInfo(edge, barcode);
        size_t frame_size = entry.GetFrameSize();
        return info.GetLeftMost() * frame_size;
    }

    /**
     *
     * @param edge
     * @param barcode
     * @return Estimated last position of the cloud defined by the barcode and the edge (not the last bin, but the last nucleotide)
     */
    size_t GetMaxPos(EdgeId edge, BarcodeId barcode) const {
        const FrameEdgeEntry<Graph> &entry = this->GetEntry(edge);
        const FrameBarcodeInfo& info = this->GetInfo(edge, barcode);
        size_t frame_size = entry.GetFrameSize();
        return info.GetRightMost() * frame_size;
    }

 private:
    using BarcodeIndexInfoExtractor<Graph, FrameEdgeEntry<Graph>>::g_;

};
typedef FrameBarcodeIndexInfoExtractorTemplate<debruijn_graph::ConjugateDeBruijnGraph> FrameBarcodeIndexInfoExtractor;
}