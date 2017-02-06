#pragma once

#include "barcode_index.hpp"

namespace barcode_index {
    class AbstractBarcodeIndexInfoExtractor {
    public:
        AbstractBarcodeIndexInfoExtractor() {}

        virtual ~AbstractBarcodeIndexInfoExtractor() {};

        virtual double GetIntersectionSizeNormalizedByUnion(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual double GetIntersectionSizeNormalizedBySecond(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual double GetIntersectionSizeNormalizedByFirst(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual size_t GetHeadBarcodeNumber(const EdgeId &edge) const = 0;

        virtual size_t GetTailBarcodeNumber(const EdgeId &edge) const = 0;

        virtual double AverageBarcodeCoverage() const = 0;

        virtual vector<int64_t> GetIntersection(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual size_t GetIntersectionSize(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual size_t GetUnionSize(const EdgeId &edge1, const EdgeId &edge2) const = 0;

        virtual bool has_barcode(const EdgeId &edge, int64_t barcode) const = 0;
    };

    template<class barcode_entry_t>
    class BarcodeIndexInfoExtractor : public AbstractBarcodeIndexInfoExtractor {
    public:
        typedef typename barcode_entry_t::barcode_distribution_t distribution_t;
        typedef typename distribution_t::key_type barcode_info_key_t;
        typedef typename distribution_t::mapped_type barcode_info_value_t;
    protected:
        shared_ptr<BarcodeIndex<barcode_entry_t>> mapper_;
        const Graph &g_;
    public:

        BarcodeIndexInfoExtractor(shared_ptr<AbstractBarcodeIndex> abstract_mapper, const Graph &g) :
                mapper_(std::dynamic_pointer_cast<BarcodeIndex<barcode_entry_t>>(abstract_mapper)),
                g_(g) {}

        double GetIntersectionSizeNormalizedByUnion(const EdgeId &edge1, const EdgeId &edge2) const override {
            if (GetUnionSize(edge1, edge2)) {
                return static_cast <double> (GetIntersectionSize(edge1, edge2)) /
                       static_cast <double> (GetUnionSize(edge1, edge2));
            }
            return 0;
        }

        virtual double GetIntersectionSizeNormalizedBySecond(const EdgeId &edge1, const EdgeId &edge2) const override {
            if (GetHeadBarcodeNumber(edge2) > 0) {
                return static_cast <double> (GetIntersectionSize(edge1, edge2)) /
                       static_cast <double> (GetHeadBarcodeNumber(edge2));
            }
            return 0;
        }

        virtual double GetIntersectionSizeNormalizedByFirst(const EdgeId &edge1, const EdgeId &edge2) const override {
            if (GetTailBarcodeNumber(edge1) > 0) {
                return static_cast <double> (GetIntersectionSize(edge1, edge2)) /
                       static_cast <double> (GetTailBarcodeNumber(edge1));
            }
            return 0;
        }


        size_t GetHeadBarcodeNumber(const EdgeId &edge) const override {
            return mapper_->GetEntryHeads(edge).Size();
        }

        size_t GetTailBarcodeNumber(const EdgeId &edge) const override {
            return mapper_->GetEntryTails(edge).Size();
        }

        double AverageBarcodeCoverage() const override {
            edge_it_helper helper(g_);
            int64_t barcodes_overall = 0;
            int64_t long_edges = 0;
            //fixme config
            size_t len_threshold = cfg::get().ts_res.edge_length_threshold;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                if (g_.length(*it) > len_threshold) {
                    long_edges++;
                    barcodes_overall += GetTailBarcodeNumber(*it);
                }
            }
            DEBUG("tails: " + std::to_string(barcodes_overall));
            DEBUG("Long edges" + long_edges);
            return static_cast <double> (barcodes_overall) / static_cast <double> (long_edges);
        }

        vector<int64_t> GetIntersection(const EdgeId &edge1, const EdgeId &edge2) const override {
            vector<int64_t> lazy_intersection;
            for (auto it = intersection_iterator_begin(edge1, edge2); it != intersection_iterator_end(edge1, edge2); ++it) {
                lazy_intersection.push_back((*it).key_);
            }
            return lazy_intersection;
        }

        size_t GetIntersectionSize(const EdgeId &edge1, const EdgeId &edge2) const override {
            auto it_tail = mapper_->GetEntryTailsIterator(edge1);
            auto it_head = mapper_->GetEntryHeadsIterator(edge2);
            return (it_tail->second).GetIntersectionSize(it_head->second);
        }

        size_t GetUnionSize(const EdgeId &edge1, const EdgeId &edge2) const override {
            auto it_tail = mapper_->GetEntryTailsIterator(edge1);
            auto it_head = mapper_->GetEntryHeadsIterator(edge2);
            return (it_tail->second).GetUnionSize(it_head->second);
        }

        bool has_barcode(const EdgeId &edge, int64_t barcode) const override {
            return mapper_->GetEntryHeads(edge).has_barcode(barcode);
        }

        typename barcode_entry_t::barcode_distribution_t::const_iterator barcode_iterator_begin(const EdgeId &edge) const {
            auto entry_it = mapper_->GetEntryHeadsIterator(edge);
            return entry_it->second.begin();
        }

        typename barcode_entry_t::barcode_distribution_t::const_iterator barcode_iterator_end(const EdgeId &edge) const {
            auto entry_it = mapper_->GetEntryHeadsIterator(edge);
            return entry_it->second.end();
        }

        const barcode_entry_t& GetEntry(const EdgeId& edge) const {
            return mapper_->GetEntryHeads(edge);
        }


        struct IntersectionData {
            const barcode_info_key_t key_;
            const barcode_info_value_t& info_first_;
            const barcode_info_value_t& info_second_;

            IntersectionData(const barcode_info_key_t key, const barcode_info_value_t& info_first, const barcode_info_value_t& info_second) :
                    key_(key), info_first_(info_first), info_second_(info_second) {}
        };

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

            //todo optimize this with lower bounds?
            const_intersection_iterator operator++() {
//                INFO(first_->first);
//                INFO(second_->first);
//                INFO(first_end_ -> first);
//                INFO(second_end_->first);
                if (first_ == first_end_ and second_ == second_end_) {
                    ++first_;
                    ++second_;
                    return *this;
                }
                if (first_->first == second_->first) {
                    if (second_ != second_end_) {
                        ++second_;
                    } else {
                        ++first_;
                    }
                }
                while (first_->first != second_->first and (first_ != first_end_ or second_ != second_end_)) {
//                    INFO("first: " << first_->first << " and");
//                    INFO("second: " << second_->first);
                    while (first_->first < second_->first and first_ != first_end_) {
                        ++first_;
//                        INFO("first: " << first_->first);
                    }
                    while (second_->first < first_->first and second_ != second_end_) {
                        ++second_;
//                        INFO("second: " << second_->first);
                    }
                    if ((first_ == first_end_ and second_->first > first_->first) or
                            (second_ == second_end_ and first_->first > second_->first)) {
                        first_ = first_end_;
                        second_ = second_end_;
                    }
                }
                if (first_->first == second_->first) {
                    return *this;
                }
                VERIFY(first_ == first_end_ and second_ == second_end_);
                if (first_->first != second_->first) {
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

            barcode_info_key_t get_first_key() {
                return first_->first;
            }
            barcode_info_key_t get_second_key() {
                return second_->first;
            }
        };

        //fixme remove end decrement
        const_intersection_iterator intersection_iterator_begin(const EdgeId& first, const EdgeId& second) const {
//            INFO(GetHeadBarcodeNumber(first));
//            INFO(GetHeadBarcodeNumber(second));
//            INFO(first.int_id());
//            INFO(second.int_id());
            if (GetHeadBarcodeNumber(first) == 0 or GetHeadBarcodeNumber(second) == 0) {
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
            if (GetHeadBarcodeNumber(first) == 0 or GetHeadBarcodeNumber(second) == 0) {
                return const_intersection_iterator(first_end, second_end, first_end, second_end);
            }
            auto first_end_copy = first_end;
            auto second_end_copy = second_end;
            const_intersection_iterator result(first_end_copy, second_end_copy, --first_end, --second_end);
            return result;
        }

    };

    class FrameBarcodeIndexInfoExtractor : public BarcodeIndexInfoExtractor<FrameEdgeEntry> {
    public:
        FrameBarcodeIndexInfoExtractor(shared_ptr<AbstractBarcodeIndex> abstract_mapper_ptr, const Graph &g) :
                BarcodeIndexInfoExtractor(abstract_mapper_ptr, g) {}

        //todo move AreEnoughSharedBarcodes to base with custom filter
        bool AreEnoughSharedBarcodes(const EdgeId &first,
                                     const EdgeId &second,
                                     size_t shared_threshold,
                                     size_t count_threshold,
                                     size_t gap_threshold) const {
            size_t current = 0;
            for (auto it = intersection_iterator_begin(first, second); it != intersection_iterator_end(first, second); ++it) {
                auto barcode = (*it).key_;
                //todo make lazy and address to info directly
                bool is_in_the_end_of_first = g_.length(first) <= gap_threshold or
                                              GetMaxPos(first, barcode) > g_.length(first) - gap_threshold;
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

        //barcode should be present on the edge
        size_t GetMinPos(const EdgeId &edge, int64_t barcode) const {
            const FrameEdgeEntry& entry = GetEntry(edge);
            const FrameBarcodeInfo& info = GetInfo(edge, barcode);
            size_t frame_size = entry.GetFrameSize();
            return info.GetLeftMost() * frame_size;
        }

        size_t GetMaxPos(const EdgeId &edge, int64_t barcode) const {
            const FrameEdgeEntry& entry = GetEntry(edge);
            const FrameBarcodeInfo& info = GetInfo(edge, barcode);
            size_t frame_size = entry.GetFrameSize();
            return info.GetRightMost() * frame_size;
        }

        const FrameBarcodeInfo& GetInfo(const EdgeId& edge, int64_t barcode) const {
            VERIFY(has_barcode(edge, barcode));
            const FrameEdgeEntry& entry = GetEntry(edge);
            return entry.get_barcode(barcode)->second;
        }

        size_t GetBarcodeLength(const EdgeId& edge, int64_t barcode) const {
            size_t max_pos = GetMaxPos(edge, barcode);
            size_t min_pos = GetMinPos(edge, barcode);
            VERIFY(max_pos >= min_pos);
            return max_pos - min_pos;
        }

        vector <size_t> GetGapDistribution(const EdgeId& edge, int64_t barcode) const {
            const FrameEdgeEntry& entry = GetEntry(edge);
            const FrameBarcodeInfo& info = GetInfo(edge, barcode);
            size_t number_of_frames = info.GetSize();
            size_t current_gap_length = 0;
            vector <size_t> result;
            for (size_t i = 0; i < number_of_frames; ++i) {
                if (not info.GetFrame(i)) {
                    ++current_gap_length;
                }
                else {
                    result.push_back(current_gap_length * entry.GetFrameSize());
                    current_gap_length = 0;
                }
            }
            std::sort(result.begin(), result.end());
            return result;
        }

//        size_t GetMedianGap(const EdgeId& edge, int64_t barcode) {
//
//
//        }
    };
}