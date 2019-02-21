#pragma once

#include <llvm/ADT/IntrusiveRefCntPtr.h>
#include "object_counter.hpp"

#include <memory>
#include <vector>

namespace pathtree {

template <typename T>
class Node : public AtomicObjectCounter<Node<T>>, public llvm::RefCountedBase<Node<T>> {
    using This = Node<T>;
    using ThisRef = llvm::IntrusiveRefCntPtr<This>;

public:
    std::vector<T> collect() const {
        std::vector<T> result;
        const This *p = this;
        while (p) {
            result.push_back(p->payload_);
            p = p->parent_.get();
        }

        return result;
    }

    Node(const T &payload, const ThisRef &parent = nullptr) : payload_{payload}, parent_{parent} {}

    static ThisRef child(const T &payload, const ThisRef &parent = nullptr) { return new This(payload, parent); }

    const auto &payload() const { return payload_; }

private:
    T payload_;
    ThisRef parent_;
};

template <class T>
using NodeRef = llvm::IntrusiveRefCntPtr<Node<T>>;

template <typename T>
NodeRef<T> make_child(const T &payload, const NodeRef<T> &parent = nullptr) {
    return Node<T>::child(payload, parent);
}

}  // namespace pathtree
