//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "common/barcode_index/barcode_index.hpp"

namespace path_extend {
namespace read_cloud {

class BarcodeScoreFunction {
  public:
    virtual ~BarcodeScoreFunction() {}

    virtual double GetScoreFromSets(const std::set<barcode_index::BarcodeId> &first,
                                    const std::set<barcode_index::BarcodeId> &second) const = 0;
};

class ContainmentIndex : BarcodeScoreFunction {
  public:
    double GetScoreFromSets(const std::set<barcode_index::BarcodeId> &first,
                            const std::set<barcode_index::BarcodeId> &second) const override {
        size_t first_size = first.size();
        size_t second_size = second.size();
        size_t min_size = std::min(first_size, second_size);
        if (min_size == 0) {
            return 0;
        }
        std::set<barcode_index::BarcodeId> intersection;
        std::set_intersection(first.begin(), first.end(), second.begin(), second.end(),
                              std::inserter(intersection, intersection.end()));
        return static_cast<double>(intersection.size()) / static_cast<double>(min_size);
    }
};
}
}