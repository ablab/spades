//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <memory>
#include <unordered_map>

namespace trie {
namespace impl {

template <typename GraphCursor>
class Node {
    using This = Node<GraphCursor>;

public:
    template <typename Collection>
    bool try_add(Collection &collection) {
        This *p = this;
        for (const GraphCursor &cursor : collection) {
            p = p->try_add_one(cursor);
            if (!p || p->is_leaf()) {
                return false;
            }
        }
        if (!p->children_.empty()) {
            return false;
        }
        p->is_leaf_ = true;
        return true;
    }

private:
    bool is_leaf_ = false;
    std::unordered_map<GraphCursor, std::unique_ptr<This>> children_;

    bool is_leaf() const { return is_leaf_; }

    This *try_add_one(const GraphCursor &cursor) {
        if (is_leaf()) {
            return nullptr;
        }

        std::unique_ptr<This> &sptr = children_[cursor];
        if (!sptr) {
            sptr.reset(new This);
        }

        return sptr.get();
    }
};

template <typename GraphCursor>
class Trie {
public:
    template <typename Collection>
    bool try_add(Collection &collection) {
        return node_.try_add(collection);
    }

private:
    Node<GraphCursor> node_;
};

}  // namespace impl

using impl::Trie;

}  // namespace trie
