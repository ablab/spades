#pragma once

#include <llvm/ADT/IntrusiveRefCntPtr.h>

#include <memory>
#include <vector>

namespace pathtrie {

namespace impl {
template <typename T>
class Node : public llvm::RefCountedBase<Node<T>> {
    using This = Node<T>;
    using ThisRef = llvm::IntrusiveRefCntPtr<This>;

public:
    // TODO Implement collect via iterators
    std::vector<T> collect() const {
        std::vector<T> result;
        const This *p = this;
        while (p) {
            result.push_back(p->data_);
            p = p->parent_.get();
        }

        return result;
    }

    static ThisRef make_child(const T &data, const ThisRef &parent = nullptr) { return new This(data, parent); }
    static ThisRef make_child(T &&data, const ThisRef &parent = nullptr) { return new This(std::move(data), parent); }
    ThisRef child(const T &data) { return make_child(data, this); }
    ThisRef child(T &&data) { return make_child(std::move(data), this); }

    auto &data() { return data_; }
    const auto &data() const { return data_; }
    bool is_root() const { return parent_ == nullptr; }

    ~Node() noexcept = default;

private:
    T data_;
    ThisRef parent_;

    Node(const T &data, const ThisRef &parent = nullptr) : data_{data}, parent_{parent} {}
    Node(T &&data, const ThisRef &parent = nullptr) : data_{std::move(data)}, parent_{parent} {}
};

template <class T>
using NodeRef = llvm::IntrusiveRefCntPtr<Node<T>>;

template <typename T>
NodeRef<T> make_root(const T &data) {
    return Node<T>::make_child(data, nullptr);
}

template <typename T>
NodeRef<T> make_root(T &&data) {
    return Node<T>::make_child(std::move(data), nullptr);
}
}

using impl::make_root;
using impl::NodeRef;
}  // namespace pathtrie
