#pragma once

#include "barcode_index.hpp"

namespace barcode_index {

    /**
     * This class provides partial interface to BarcodeIndexInfoExtractor.
     */
    class AbstractBarcodeIndexInfoExtractor {
    public:
        AbstractBarcodeIndexInfoExtractor() {}

        virtual ~AbstractBarcodeIndexInfoExtractor() {};

        virtual size_t GetNumberOfBarcodes(const EdgeId &edge) const = 0;

        virtual vector<BarcodeId> GetSharedBarcodes(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual size_t GetNumberOfSharedBarcodes(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual bool HasBarcode(const EdgeId &edge, const BarcodeId& barcode) const = 0;

        virtual double AverageBarcodeCoverage() const = 0;

        virtual double GetIntersectionSizeNormalizedByUnion(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual double GetIntersectionSizeNormalizedBySecond(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual size_t GetUnionSize(const EdgeId &edge1, const EdgeId &edge2) const = 0;
    };

    /**
     * BarcodeIndexInfoExtractor extracts useful statistics from BarcodeIndex
     */
    template<class barcode_entry_t>
    class BarcodeIndexInfoExtractor : public AbstractBarcodeIndexInfoExtractor {
    public:
        typedef typename barcode_entry_t::barcode_distribution_t distribution_t;
        typedef typename barcode_entry_t::barcode_info_t barcode_info_t;
        typedef typename distribution_t::key_type barcode_info_key_t;
        typedef typename distribution_t::mapped_type barcode_info_value_t;
        typedef typename distribution_t::value_type barcode_info_pair_t;
    protected:
        shared_ptr<barcode_index::BarcodeIndex<barcode_entry_t>> mapper_;
        const Graph &g_;
    public:
        BarcodeIndexInfoExtractor(shared_ptr<AbstractBarcodeIndex> abstract_mapper, const Graph &g) :
                mapper_(std::dynamic_pointer_cast<barcode_index::BarcodeIndex<barcode_entry_t>>(abstract_mapper)),
                g_(g) {}

        /**
         *
         * @param edge
         * @return Number of barcodes contained by the edge
         */
        size_t GetNumberOfBarcodes(const EdgeId &edge) const override {
            return mapper_->GetEntry(edge).Size();
        }

        /**
         *
         * @param edge1
         * @param edge2
         * @return List of barcodes shared by edge1 and edge2
         */
        vector<BarcodeId> GetSharedBarcodes (const EdgeId &edge1, const EdgeId &edge2) const override {
            vector<BarcodeId> intersection;
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
        size_t GetNumberOfSharedBarcodes (const EdgeId &edge1, const EdgeId &edge2) const override {
            auto it_tail = mapper_->GetEntryTailsIterator(edge1);
            auto it_head = mapper_->GetEntryHeadsIterator(edge2);
            return (it_tail->second).GetIntersectionSize(it_head->second);
        }

        /**
         * @param edge
         * @param barcode
         * @return True if the edge contains the barcode
         */
        bool HasBarcode(const EdgeId &edge, const BarcodeId& barcode) const override {
            return mapper_->GetEntry(edge).has_barcode(barcode);
        }

        /**
         *
         * @return Average number of barcodes contained on the long edges of the graph
         */
        double AverageBarcodeCoverage() const override {
            edge_it_helper helper(g_);
            size_t barcodes_overall = 0;
            size_t long_edges = 0;
            size_t len_threshold = cfg::get().ts_res.edge_length_threshold;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                if (g_.length(*it) > len_threshold) {
                    long_edges++;
                    barcodes_overall += GetNumberOfBarcodes(*it);
                }
            }
            DEBUG("tails: " + std::to_string(barcodes_overall));
            DEBUG("Long edges" + long_edges);
            return static_cast <double> (barcodes_overall) / static_cast <double> (long_edges);
        }

        double GetIntersectionSizeNormalizedByUnion(const EdgeId &edge1, const EdgeId &edge2) const override {
            if (GetUnionSize(edge1, edge2)) {
                return static_cast <double> (GetNumberOfSharedBarcodes(edge1, edge2)) /
                       static_cast <double> (GetUnionSize(edge1, edge2));
            }
            return 0;
        }

        virtual double GetIntersectionSizeNormalizedBySecond(const EdgeId &edge1, const EdgeId &edge2) const override {
            if (GetNumberOfBarcodes(edge2) > 0) {
                return static_cast <double> (GetNumberOfSharedBarcodes(edge1, edge2)) /
                       static_cast <double> (GetNumberOfBarcodes(edge2));
            }
            return 0;
        }

        size_t GetUnionSize(const EdgeId &edge1, const EdgeId &edge2) const override {
            auto it_tail = mapper_->GetEntryTailsIterator(edge1);
            auto it_head = mapper_->GetEntryHeadsIterator(edge2);
            return (it_tail->second).GetUnionSize(it_head->second);
        }

        vector<BarcodeId> GetBarcodes(const EdgeId& edge) const {
            vector <BarcodeId> result;
            auto copy_barcode_id = [&result](const barcode_info_pair_t& entry)
                        {result.push_back(entry.first); };
            std::for_each(barcode_iterator_begin(edge), barcode_iterator_end(edge), copy_barcode_id);
            return result;
        }


        typename distribution_t::const_iterator barcode_iterator_begin(const EdgeId &edge) const {
            auto entry_it = mapper_->GetEntryHeadsIterator(edge);
            return entry_it->second.begin();
        }

        typename distribution_t::const_iterator barcode_iterator_end(const EdgeId &edge) const {
            auto entry_it = mapper_->GetEntryHeadsIterator(edge);
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

            IntersectionData(const barcode_info_key_t key, const barcode_info_value_t& info_first, const barcode_info_value_t& info_second) :
                    key_(key), info_first_(info_first), info_second_(info_second) {}
        };

        /**
         * Iterator over shared barcodes of two edges.
         * Dereferencing returns proxy object of type IntersectionData
         * @note Since it is not an iterator over container there is no -> operator.
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

            //todo optimize with lower bounds
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
                    DEBUG("first: " << get_first_key() << " and");
                    DEBUG("second: " << get_second_key());
                    while (get_first_key() < get_second_key() and first_ != first_end_) {
                        ++first_;
                        DEBUG("first: " << get_first_key());
                    }
                    while (get_second_key() < get_first_key() and second_ != second_end_) {
                        ++second_;
                        DEBUG("second: " << get_second_key());
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
                return not(*this == other);
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
        const_intersection_iterator intersection_iterator_begin(const EdgeId& first, const EdgeId& second) const {
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

        const_intersection_iterator intersection_iterator_end(const EdgeId& first, const EdgeId& second) const {
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
        const barcode_info_t& GetInfo(const EdgeId& edge, const BarcodeId& barcode) const {
            VERIFY(HasBarcode(edge, barcode));
            const barcode_entry_t& entry = GetEntry(edge);
            return entry.get_barcode(barcode)->second;
        }

        const barcode_entry_t& GetEntry(const EdgeId& edge) const {
            return mapper_->GetEntry(edge);
        }

    };

    /**
     * Specialization of BarcodeIndexInfoExtractor for FrameBarcodeInfo type.
     * @see FrameBarcodeInfo
     */
    class FrameBarcodeIndexInfoExtractor : public BarcodeIndexInfoExtractor<FrameEdgeEntry> {
    public:
        FrameBarcodeIndexInfoExtractor(shared_ptr<AbstractBarcodeIndex> abstract_mapper_ptr, const Graph &g) :
                BarcodeIndexInfoExtractor(abstract_mapper_ptr, g) {}

        /**
         *
         * @param edge
         * @param barcode
         * @return number of barcoded reads aligned to the edge
         */
        size_t GetNumberOfReads(const EdgeId &edge, const BarcodeId &barcode) const {
            return GetInfo(edge, barcode).GetCount();
        }

        /**
         *
         * @param edge
         * @param barcode
         * @return leftmost barcoded bin of the edge
         */
        size_t GetLeftBin(const EdgeId &edge, const BarcodeId &barcode) const {
            return GetInfo(edge, barcode).GetLeftMost();
        }

        /**
         *
         * @param edge
         * @param barcode
         * @return rightmost barcoded bin of the edge
         */
        size_t GetRightBin(const EdgeId &edge, const BarcodeId &barcode) const {
            return GetInfo(edge, barcode).GetRightMost();
        }


        /**
         *
         * @param edge
         * @param barcode
         * @return bitset representing barcoded bins of the edge
         */
        const boost::dynamic_bitset<>& GetBitSet(const EdgeId& edge, const BarcodeId& barcode) const {
            return GetInfo(edge, barcode).GetBitSet();
        }

        /**
         * @param edge
         * @return length of the bin
         */
        size_t GetBinLength(const EdgeId &edge) const {
            return GetEntry(edge).GetFrameSize();
        }

        /**
         *
         * @param edge
         * @return number of bins on the edge
         */
        size_t GetNumberOfBins(const EdgeId& edge) const {
            return GetEntry(edge).GetNumberOfFrames();
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
        bool AreEnoughSharedBarcodesWithFilter (const EdgeId &first,
                                                const EdgeId &second,
                                                size_t shared_threshold,
                                                size_t count_threshold,
                                                size_t gap_threshold) const {
            size_t current = 0;
            for (auto it = intersection_iterator_begin(first, second); it != intersection_iterator_end(first, second); ++it) {
                BarcodeId barcode = (*it).key_;
                //todo make lazy and address to info directly
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
        size_t CountSharedBarcodesWithFilter (const EdgeId &first,
                                             const EdgeId &second,
                                             size_t count_threshold,
                                             size_t gap_threshold) const {
            size_t current = 0;
            for (auto it = intersection_iterator_begin(first, second); it != intersection_iterator_end(first, second); ++it) {
                auto barcode = (*it).key_;
                bool is_in_the_end_of_first = g_.length(first) <= gap_threshold or
                                                     GetMaxPos(first, barcode) + gap_threshold > g_.length(first);
                bool is_in_the_beginning_of_second = g_.length(second) <= gap_threshold or
                                                     GetMinPos(second, barcode) < gap_threshold;
                bool enough_count = (*it).info_first_.GetCount() >= count_threshold and
                                    (*it).info_second_.GetCount() >= count_threshold;
                if (is_in_the_end_of_first and is_in_the_beginning_of_second and enough_count) {
                    ++current;
                }
            }
            return current;
        }

        /**
         * @param first first edge
         * @param second second edge
         * @param count_threshold edge contains barcode iff there are at least count_threshold reads aligned to the edge
         * @param gap_threshold clouds located at the beginning of the first or at the end of the second edge are discarded.
         *      Cloud is located in the beginning of the edge if it is not aligned to the last gap_threshold nucleotides of the edge.
         * @return list of barcodes which pass requirements determined by count_threshold and gap_threshold.
         */
        vector<BarcodeId> GetSharedBarcodesWithFilter (const EdgeId &first,
                                              const EdgeId &second,
                                              size_t count_threshold,
                                              size_t gap_threshold) const {
            vector<BarcodeId> result;
            for (auto it = intersection_iterator_begin(first, second); it != intersection_iterator_end(first, second); ++it) {
                auto barcode = (*it).key_;
                bool is_in_the_end_of_first = g_.length(first) <= gap_threshold or
                                                     GetMaxPos(first, barcode) + gap_threshold > g_.length(first);
                bool is_in_the_beginning_of_second = g_.length(second) <= gap_threshold or
                                                     GetMinPos(second, barcode) < gap_threshold;
                bool enough_count = (*it).info_first_.GetCount() >= count_threshold and
                                    (*it).info_second_.GetCount() >= count_threshold;
                if (is_in_the_end_of_first and is_in_the_beginning_of_second and enough_count) {
                    result.push_back(barcode);
                }
            }
            return result;
        }

        /**
         *
         * @param edge
         * @param barcode
         * @return Estimated first position of the cloud defined by the barcode and the edge (not the first bin, but the first nucleotide)
         */
        size_t GetMinPos(const EdgeId &edge, const BarcodeId& barcode) const {
            const FrameEdgeEntry& entry = GetEntry(edge);
            const FrameBarcodeInfo& info = GetInfo(edge, barcode);
            size_t frame_size = entry.GetFrameSize();
            return info.GetLeftMost() * frame_size;
        }

        /**
         *
         * @param edge
         * @param barcode
         * @return Estimated last position of the cloud defined by the barcode and the edge (not the last bin, but the last nucleotide)
         */
        size_t GetMaxPos(const EdgeId &edge, const BarcodeId& barcode) const {
            const FrameEdgeEntry& entry = GetEntry(edge);
            const FrameBarcodeInfo& info = GetInfo(edge, barcode);
            size_t frame_size = entry.GetFrameSize();
            return info.GetRightMost() * frame_size;
        }


        size_t GetBarcodeLength(const EdgeId& edge, const BarcodeId& barcode) const {
            size_t max_pos = GetMaxPos(edge, barcode);
            size_t min_pos = GetMinPos(edge, barcode);
            VERIFY(max_pos >= min_pos);
            return max_pos - min_pos;
        }
    };
}