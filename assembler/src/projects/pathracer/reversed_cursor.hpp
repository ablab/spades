//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <vector>

template <class GraphCursor>
class ReversedGraphCursor : public GraphCursor {
public:
    using Context = typename GraphCursor::Context;
    using GraphCursor::GraphCursor;
    using GraphCursor::operator=;
    // "an inherited constructor is not a candidate for initialization from an expression of the same or derived type"
    ReversedGraphCursor(const GraphCursor& other) : GraphCursor(other) {}
    ReversedGraphCursor(GraphCursor&& other) : GraphCursor(std::move(other)) {}

    ReversedGraphCursor() = default;
    ReversedGraphCursor(const ReversedGraphCursor&) = default;
    ReversedGraphCursor(ReversedGraphCursor&&) = default;
    ReversedGraphCursor& operator=(const ReversedGraphCursor&) = default;
    ReversedGraphCursor& operator=(ReversedGraphCursor&&) = default;

    std::vector<ReversedGraphCursor> next(typename GraphCursor::Context context) const {
        auto result = GraphCursor::prev(context);
        return std::vector<ReversedGraphCursor>(std::cbegin(result), std::cend(result));
    }

    std::vector<ReversedGraphCursor> prev(typename GraphCursor::Context context) const {
        auto result = GraphCursor::next(context);
        return std::vector<ReversedGraphCursor>(std::cbegin(result), std::cend(result));
    }
};

namespace std {
template <typename GraphCursor>
struct hash<ReversedGraphCursor<GraphCursor>> : public hash<GraphCursor> {};
}  // namespace std
