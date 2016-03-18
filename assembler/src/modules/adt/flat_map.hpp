#ifndef __ADT_FLAT_MAP_HPP__
#define __ADT_FLAT_MAP_HPP__

#pragma once

#include <vector>
#include <algorithm>
#include <functional>

namespace adt {

template<typename K, typename V, typename Comp = std::less<K>, typename Allocator = std::allocator<std::pair<K, V> > >
struct flat_map {
    typedef K key_type;
    typedef V mapped_type;
    typedef std::pair<K, V> value_type;
    typedef Comp key_compare;
    struct value_compare : std::binary_function<value_type, value_type, bool> {
        bool operator()(const value_type & lhs, const value_type & rhs) const {
            return key_compare()(lhs.first, rhs.first);
        }
    };
    typedef Allocator allocator_type;
    typedef V& reference;
    typedef const V& const_reference;
    typedef typename std::allocator_traits<allocator_type>::pointer pointer;
    typedef typename std::allocator_traits<allocator_type>::const_pointer const_pointer;
    typedef std::vector<value_type, allocator_type> container_type;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;
    typedef typename container_type::reverse_iterator reverse_iterator;
    typedef typename container_type::const_reverse_iterator const_reverse_iterator;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::size_type size_type;

    flat_map() = default;
    template<typename It>
    flat_map(It begin, It end) { insert(begin, end); }
    flat_map(std::initializer_list<value_type> init)
            : flat_map(init.begin(), init.end()) {}

    iterator                begin()              {    return data_.begin();    }
    iterator                end()                {    return data_.end();      }
    const_iterator          begin()     const    {    return data_.begin();    }
    const_iterator          end()       const    {    return data_.end();      }
    const_iterator          cbegin()    const    {    return data_.cbegin();   }
    const_iterator          cend()      const    {    return data_.cend();     }
    reverse_iterator        rbegin()             {    return data_.rbegin();   }
    reverse_iterator        rend()               {    return data_.rend();     }
    const_reverse_iterator  rbegin()    const    {    return data_.rbegin();   }
    const_reverse_iterator  rend()      const    {    return data_.rend();     }
    const_reverse_iterator  crbegin()   const    {    return data_.crbegin();  }
    const_reverse_iterator  crend()     const    {    return data_.crend();    }

    bool empty() const { return data_.empty(); }
    size_type size() const { return data_.size(); }
    size_type max_size() const { return data_.max_size(); }
    size_type capacity() const { return data_.capacity(); }
    void reserve(size_type size) {data_.reserve(size); }
    void shrink_to_fit() { data_.shrink_to_fit(); }
    size_type bytes_used() const { return capacity() * sizeof(value_type) + sizeof(data_); }

    mapped_type & operator[](const key_type &key) {
        KeyOrValueCompare comp;
        auto lower = lower_bound(key);
        if (lower == end() || comp(key, *lower))
            return data_.emplace(lower, key, mapped_type())->second;
        else
            return lower->second;
    }
    mapped_type & operator[](key_type &&key) {
        KeyOrValueCompare comp;
        auto lower = lower_bound(key);
        if (lower == end() || comp(key, *lower))
            return data_.emplace(lower, std::move(key), mapped_type())->second;
        else
            return lower->second;
    }

    std::pair<iterator, bool> insert(value_type &&value) {
        return emplace(std::move(value));
    }
    std::pair<iterator, bool> insert(const value_type &value) {
        return emplace(value);
    }
    iterator insert(const_iterator hint, value_type &&value) {
        return emplace_hint(hint, std::move(value));
    }
    iterator insert(const_iterator hint, const value_type &value) {
        return emplace_hint(hint, value);
    }

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

        // Insert the remaining elements that didn't fit by calling this function recursively.
        return insert(begin, end);
    }
    void insert(std::initializer_list<value_type> il) {
        insert(il.begin(), il.end());
    }
    iterator erase(iterator it) {
        return data_.erase(it);
    }
    iterator erase(const_iterator it) {
        return erase(iterator_const_cast(it));
    }
    size_type erase(const key_type &key) {
        auto found = find(key);
        if (found == end())
            return 0;
        erase(found);
        return 1;
    }
    iterator erase(const_iterator first, const_iterator last) {
        return data_.erase(iterator_const_cast(first), iterator_const_cast(last));
    }
    void swap(flat_map & other) {
        data_.swap(other.data);
    }
    void clear() {
        data_.clear();
    }
    template<typename First, typename... Args>
    std::pair<iterator, bool> emplace(First &&first, Args &&... args) {
        KeyOrValueCompare comp;
        auto lower_bound = std::lower_bound(data_.begin(), data_.end(), first, comp);
        if (lower_bound == data_.end() || comp(first, *lower_bound))
            return { data_.emplace(lower_bound, std::forward<First>(first), std::forward<Args>(args)...), true };
        else
            return { lower_bound, false };
    }
    std::pair<iterator, bool> emplace() {
        return emplace(value_type());
    }
    template<typename First, typename... Args>
    iterator emplace_hint(const_iterator hint, First &&first, Args &&... args) {
        KeyOrValueCompare comp;
        if (hint == cend() || comp(first, *hint)) {
            if (hint == cbegin() || comp(*(hint - 1), first))
                return data_.emplace(iterator_const_cast(hint), std::forward<First>(first), std::forward<Args>(args)...);
            else
                return emplace(std::forward<First>(first), std::forward<Args>(args)...).first;
        } else if (!comp(*hint, first)) {
            return begin() + (hint - cbegin());
        } else {
            return emplace(std::forward<First>(first), std::forward<Args>(args)...).first;
        }
    }
    iterator emplace_hint(const_iterator hint) {
        return emplace_hint(hint, value_type());
    }

    key_compare key_comp() const {
        return key_compare();
    }
    value_compare value_comp() const {
        return value_compare();
    }

    template<typename T>
    iterator find(const T &key) {
        return binary_find(begin(), end(), key, KeyOrValueCompare());
    }
    template<typename T>
    const_iterator find(const T &key) const {
        return binary_find(begin(), end(), key, KeyOrValueCompare());
    }
    template<typename T>
    size_type count(const T &key) const {
        return std::binary_search(begin(), end(), key, KeyOrValueCompare()) ? 1 : 0;
    }
    template<typename T>
    iterator lower_bound(const T &key) {
        return std::lower_bound(begin(), end(), key, KeyOrValueCompare());
    }
    template<typename T>
    const_iterator lower_bound(const T & key) const {
        return std::lower_bound(begin(), end(), key, KeyOrValueCompare());
    }
    template<typename T>
    iterator upper_bound(const T & key) {
        return std::upper_bound(begin(), end(), key, KeyOrValueCompare());
    }
    template<typename T>
    const_iterator upper_bound(const T &key) const {
        return std::upper_bound(begin(), end(), key, KeyOrValueCompare());
    }
    template<typename T>
    std::pair<iterator, iterator> equal_range(const T &key) {
        return std::equal_range(begin(), end(), key, KeyOrValueCompare());
    }
    template<typename T>
    std::pair<const_iterator, const_iterator> equal_range(const T &key) const {
        return std::equal_range(begin(), end(), key, KeyOrValueCompare());
    }
    allocator_type get_allocator() const {
        return data_.get_allocator();
    }

    bool operator==(const flat_map &other) const {
        return data_ == other.data_;
    }
    bool operator!=(const flat_map &other) const {
        return !(*this == other);
    }
    bool operator<(const flat_map &other) const {
        return data_ < other.data_;
    }
    bool operator>(const flat_map &other) const {
        return other < *this;
    }
    bool operator<=(const flat_map &other) const {
        return !(other < *this);
    }
    bool operator>=(const flat_map &other) const {
        return !(*this < other);
    }

  private:
    container_type data_;

    iterator iterator_const_cast(const_iterator it) {
        return begin() + (it - cbegin());
    }

    struct KeyOrValueCompare {
        bool operator()(const key_type &lhs, const key_type &rhs) const {
            return key_compare()(lhs, rhs);
        }
        bool operator()(const key_type &lhs, const value_type &rhs) const {
            return key_compare()(lhs, rhs.first);
        }
        template<typename T>
        bool operator()(const key_type &lhs, const T &rhs) const {
            return key_compare()(lhs, rhs);
        }
        template<typename T>
        bool operator()(const T &lhs, const key_type &rhs) const {
            return key_compare()(lhs, rhs);
        }
        bool operator()(const value_type &lhs, const key_type &rhs) const {
            return key_compare()(lhs.first, rhs);
        }
        bool operator()(const value_type &lhs, const value_type &rhs) const {
            return key_compare()(lhs.first, rhs.first);
        }
        template<typename T>
        bool operator()(const value_type &lhs, const T &rhs) const {
            return key_compare()(lhs.first, rhs);
        }
        template<typename T>
        bool operator()(const T &lhs, const value_type &rhs) const {
            return key_compare()(lhs, rhs.first);
        }
    };

    // like std::binary_search, but returns the iterator to the element
    // if it was found, and returns end otherwise
    template<typename It, typename T, typename Compare>
    static It binary_find(It begin, It end, const T & value, const Compare & cmp) {
        auto lower_bound = std::lower_bound(begin, end, value, cmp);
        if (lower_bound == end || cmp(value, *lower_bound))
            return end;
        else
            return lower_bound;
    }
};

template<typename K, typename V, typename C, typename A>
void swap(flat_map<K, V, C, A> & lhs, flat_map<K, V, C, A> & rhs) {
    lhs.swap(rhs);
}

}

#endif
