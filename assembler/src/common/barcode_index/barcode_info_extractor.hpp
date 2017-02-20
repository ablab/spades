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
    protected:
        typedef typename barcode_entry_t::barcode_distribution_t distribution_t;
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
            auto it_tail = mapper_->GetEntryTailsIterator(edge1);
            auto it_head = mapper_->GetEntryHeadsIterator(edge2);
            return (it_tail->second).GetIntersection(it_head->second);
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
//
//        class const_intersection_iterator : public std::iterator<std::input_iterator_tag, distribution_t::value_type> {
//        public:
//            const_intersection_iterator(distribution_t::const_iterator first, distribution_t::const_iterator second) :
//                    first_(first), second_(second) {}
//
//            const_intersection_iterator operator++() {
//
//            }
//
//        private:
//            distribution_t::const_iterator first_;
//            distribution_t::const_iterator second_;
//        };


    };

    class FrameBarcodeIndexInfoExtractor : public BarcodeIndexInfoExtractor<FrameEdgeEntry> {
    public:
        FrameBarcodeIndexInfoExtractor(shared_ptr<AbstractBarcodeIndex> abstract_mapper_ptr, const Graph &g) :
                BarcodeIndexInfoExtractor(abstract_mapper_ptr, g) {}

        size_t GetNumberOfSharedWithFilter(const EdgeId &first, const EdgeId &second, size_t gap_threshold) const {
            //fixme implement intersection iterator
            auto barcodes = GetIntersection(first, second);
            size_t result = 0;
            for (auto barcode: barcodes) {
                if (g_.length(first) <= gap_threshold or
                        GetMaxPos(first, barcode) > g_.length(first) - gap_threshold) {
                    if (g_.length(second) <= gap_threshold or
                            GetMinPos(second, barcode) < gap_threshold) {
                        ++result;
                    }
                }
            }
            return result;
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