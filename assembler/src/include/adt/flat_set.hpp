#ifndef __ADT_FLAT_SET_HPP__
#define __ADT_FLAT_SET_HPP__

#pragma once

#include <vector>
#include <algorithm>
#include <type_traits>
#include <functional>

namespace adt {

template<typename T, typename Comp = std::less<T>, template<typename, typename...> class Container = std::vector >
struct flat_set {
    typedef T key_type;
    typedef T value_type;
    typedef Comp key_compare;
    typedef Comp value_compare;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef Container<value_type> container_type;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;
    typedef typename container_type::reverse_iterator reverse_iterator;
    typedef typename container_type::const_reverse_iterator const_reverse_iterator;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::size_type size_type;

    flat_set() = default;
    template<typename It>
    flat_set(It begin, It end) {
        insert(begin, end);
    }
    flat_set(std::initializer_list<value_type> init)
            : flat_set(init.begin(), init.end()) { }

    iterator                begin()         { return data_.begin();   }
    iterator                end()           { return data_.end();     }
    const_iterator          begin()   const { return data_.begin();   }
    const_iterator          end()     const { return data_.end();     }
    const_iterator          cbegin()  const { return data_.cbegin();  }
    const_iterator          cend()    const { return data_.cend();    }
    reverse_iterator        rbegin()        { return data_.rbegin();  }
    reverse_iterator        rend()          { return data_.rend();    }
    const_reverse_iterator  rbegin()  const { return data_.rbegin();  }
    const_reverse_iterator  rend()    const { return data_.rend();    }
    const_reverse_iterator  crbegin() const { return data_.crbegin(); }
    const_reverse_iterator  crend()   const { return data_.crend();   }

    bool empty() const { return data_.empty(); }
    size_type size() const { return data_.size(); }
    size_type max_size() const { return data_.max_size(); }
    size_type capacity() const { return data_.capacity(); }
    void reserve(size_type size) { data_.reserve(size); }
    void shrink_to_fit() { data_.shrink_to_fit(); }
    size_type bytes_used() const { return capacity() * sizeof(value_type) + sizeof(data_); }

    std::pair<iterator, bool> insert(value_type && value) { return emplace(std::move(value)); }
    std::pair<iterator, bool> insert(const value_type & value) { return emplace(value); }
    iterator insert(const_iterator hint, value_type && value) { return emplace_hint(hint, std::move(value)); }
    iterator insert(const_iterator hint, const value_type & value) { return emplace_hint(hint, value); }
    void insert(std::initializer_list<value_type> il) { insert(il.begin(), il.end()); }

    template<typename It>
    void insert(It begin, It end) {
        // If we need to increase the capacity, utilize this fact and emplace
        // the stuff.
        for (; begin != end && size() == capacity(); ++begin) {
            emplace(*begin);
        }
        if (begin == end)
            return;
        // If we don't need to increase capacity, then we can use a more efficient
        // insert method where everything is just put in the same vector
        // and then merge in place.
        size_type size_before = data_.size();
        try {
            for (size_t i = capacity(); i > size_before && begin != end; --i, ++begin) {
                data_.emplace_back(*begin);
            }
        } catch(...) {
            // If emplace_back throws an exception, the easiest way to make sure
            // that our invariants are still in place is to resize to the state
            // we were in before
            for (size_t i = data_.size(); i > size_before; --i) {
                data_.pop_back();
            }
            throw;
        }
        value_compare comp;
        auto mid = data_.begin() + size_before;
        std::stable_sort(mid, data_.end(), comp);
        std::inplace_merge(data_.begin(), mid, data_.end(), comp);
        data_.erase(std::unique(data_.begin(), data_.end(), std::not2(comp)), data_.end());
        // Make sure that we inserted at least one element before
        // recursing. Otherwise we'd recurse too often if we were to insert the
        // same element many times
        if (data_.size() == size_before) {
            for (; begin != end; ++begin) {
                if (emplace(*begin).second) {
                    ++begin;
                    break;
                }
            }
        }

        // insert the remaining elements that didn't fit by calling this function recursively
        // this will recurse log(n) times where n is std::distance(begin, end)
        return insert(begin, end);
    }
    iterator erase(iterator it) { return data_.erase(it); }
    iterator erase(const_iterator it) { return erase(iterator_const_cast(it)); }
    size_type erase(const value_type &val) {
        auto found = find(val);
        if (found == end()) return 0;
        erase(found);
        return 1;
    }
    iterator erase(const_iterator first, const_iterator last) {
        return data_.erase(iterator_const_cast(first), iterator_const_cast(last));
    }

    void swap(flat_set & other) { data_.swap(other.data); }
    void clear() { data_.clear(); }

    template<typename First, typename... Args>
    std::pair<iterator, bool> emplace(First && first, Args &&... args) {
        Comp comp;
        auto lower_bound = std::lower_bound(data_.begin(), data_.end(), first, comp);
        if (lower_bound == data_.end() || comp(first, *lower_bound))
            return { data_.emplace(lower_bound, std::forward<First>(first), std::forward<Args>(args)...), true };
        else
            return { lower_bound, false };
    }
    std::pair<iterator, bool> emplace() { return emplace(value_type()); }
    template<typename First, typename... Args>
    iterator emplace_hint(const_iterator hint, First && first, Args &&... args) {
        Comp comp;
        if (hint == cend() || comp(first, *hint)) {
            if (hint == cbegin() || comp(*(hint - 1), first))
                return data_.emplace(iterator_const_cast(hint), std::forward<First>(first), std::forward<Args>(args)...);
            else
                return emplace(std::forward<First>(first), std::forward<Args>(args)...).first;
        } else if (!comp(*hint, first)) {
            return begin() + (hint - cbegin());
        }

        return emplace(std::forward<First>(first), std::forward<Args>(args)...).first;
    }
    iterator emplace_hint(const_iterator hint) { return emplace_hint(hint, value_type()); }

    key_compare key_comp() const { return key_compare(); }
    value_compare value_comp() const { return value_compare(); }

    iterator find(const value_type &key) {
        return binary_find(begin(), end(), key, Comp());
    }
    const_iterator find(const value_type &key) const {
        return binary_find(begin(), end(), key, Comp());
    }
    size_type count(const value_type &key) const {
        return std::binary_search(begin(), end(), key, Comp()) ? 1 : 0;
    }
    iterator lower_bound(const value_type &key) {
        return std::lower_bound(begin(), end(), key, Comp());
    }
    const_iterator lower_bound(const value_type &key) const {
        return std::lower_bound(begin(), end(), key, Comp());
    }
    iterator upper_bound(const value_type &key) {
        return std::upper_bound(begin(), end(), key, Comp());
    }
    const_iterator upper_bound(const value_type &key) const {
        return std::upper_bound(begin(), end(), key, Comp());
    }
    std::pair<iterator, iterator> equal_range(const value_type &key) {
        return std::equal_range(begin(), end(), key, Comp());
    }
    std::pair<const_iterator, const_iterator> equal_range(const value_type &key) const {
        return std::equal_range(begin(), end(), key, Comp());
    }

    bool operator==(const flat_set &other) const {
        return data_ == other.data_;
    }
    bool operator!=(const flat_set &other) const {
        return !(*this == other);
    }
    bool operator<(const flat_set &other) const {
        return data_ < other.data_;
    }
    bool operator>(const flat_set &other) const {
        return other < *this;
    }
    bool operator<=(const flat_set &other) const {
        return !(other < *this);
    }
    bool operator>=(const flat_set &other) const {
        return !(*this < other);
    }

  private:
    container_type data_;

    iterator iterator_const_cast(const_iterator it) {
        return begin() + (it - cbegin());
    }

    // like std::binary_search, but returns the iterator to the element
    // if it was found, and returns end otherwise
    template<typename It, typename Compare>
    static It binary_find(It begin, It end, const value_type &value, const Compare &cmp) {
        auto lower_bound = std::lower_bound(begin, end, value, cmp);
        if (lower_bound == end || cmp(value, *lower_bound))
            return end;
        else
            return lower_bound;
    }
};

template<typename V, typename C, template<typename, typename...> class Container>
void swap(flat_set<V, C, Container> & lhs, flat_set<V, C, Container> & rhs) {
    lhs.swap(rhs);
}

}

#endif // __ADT_FLAT_SET_HPP__
