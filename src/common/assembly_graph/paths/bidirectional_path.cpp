//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/paths/bidirectional_path.hpp"
#include <vector>
#include <string>

namespace path_extend {

std::atomic<uint64_t> BidirectionalPath::path_id_{0};

std::vector<std::string> BidirectionalPath::PrintLines() const {
    auto as_str = str();
    utils::trim(as_str);
    
    auto res = utils::split(as_str, "\n", true);

    return { res.begin(), res.end() };
}

void BidirectionalPath::PrintDEBUG() const {
    DEBUG_EXPR( // This will not execute PrintLines() if DEBUG level is not reached.
        for (const auto& s: PrintLines()) {
            DEBUG(s);
        }
    );
}

void BidirectionalPath::PrintINFO() const {
    for (const auto& s: PrintLines()) {
        INFO(s);
    }
}

void BidirectionalPath::Print(std::ostream &os) const {
    if (Empty())
        return;

    os << "Path " << GetId() << "\n";
    os << "Length " << Length() << "\n";
    os << "Weight " << weight_ << "\n";
    os << "#, edge (length), gap info, total length, total length from start" << "\n";
    for (size_t i = 0; i < Size(); ++i) {
        os << i << ", " << g_.str(At(i))
           << ", " << GapAt(i)
           << ", " << LengthAt(i)
           << ", " << ((Length() < LengthAt(i)) ? 0 : Length() - LengthAt(i)) << "\n";
    }
}

std::string BidirectionalPath::str() const {
    std::stringstream ss;
    Print(ss);
    return ss.str();
}


}
