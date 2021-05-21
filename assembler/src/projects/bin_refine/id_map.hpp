//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "adt/iterator_range.hpp"
#include <boost/iterator/iterator_facade.hpp>

#include <functional>
#include <type_traits>
#include <vector>

namespace adt {
template<typename V, typename K = uint64_t>
class id_map {
  public:
    typedef K key_type;
    typedef V mapped_type;

    typedef std::vector<mapped_type> container_type;

    template<bool IsConst = false>
    class iterator : public boost::iterator_facade<iterator<IsConst>,
                                                   std::conditional_t<IsConst, const mapped_type, mapped_type>,
                                                   boost::forward_traversal_tag> {
        typedef typename std::conditional_t<IsConst, const mapped_type, mapped_type> value_type;
        typedef typename std::conditional_t<IsConst, const mapped_type&, mapped_type&> reference_type;
        typedef typename std::conditional_t<IsConst, const container_type, container_type> data_type;

      public:
        static constexpr size_t NPOS = -1ULL;

        iterator(size_t start,
                 const std::vector<bool> &map,
                 data_type &data)
                : map_(map), data_(data), cur_(start) {
            if (cur_ != NPOS && !map_.get()[cur_])
                cur_ = next_occupied(cur_);
        }

        key_type key() const {
            return cur_;
        }

        reference_type value() const {
            return dereference();
        }
        
      private:
        friend class boost::iterator_core_access;

        reference_type& dereference() const {
            return data_.get()[cur_];
        }

        size_t next_occupied(uint64_t n) const {
            for (size_t i = n + 1; i < map_.get().size(); ++i) {
                if (map_.get()[i])
                    return i;
            }
            return NPOS;
        }

        void increment() {
            if (cur_ == NPOS)
                return;

            cur_ = next_occupied(cur_);
        }

        bool equal(const iterator &other) const {
            return cur_ == other.cur_;
        }

      private:
        std::reference_wrapper<const std::vector<bool>> map_;
        std::reference_wrapper<data_type> data_;
        size_t cur_;
    };

    typedef iterator<true> const_iterator;

    id_map(size_t max_id) {
        data_.resize(max_id);
        occ_map_.resize(max_id, false);
    }
    
    mapped_type& at(const key_type key) {
        return data_.at(key.int_id());
    }
    const mapped_type& at(const key_type key) const {
        return data_.at(key.int_id());
    }

    mapped_type& operator[](const key_type key) noexcept {
        size_t id = key.int_id();
        if (!occ_map_[id]) occ_map_[id] = true;
        return data_[id];
    }

    const mapped_type& operator[](const key_type key) const noexcept {
        return data_[key.int_id()];
    }

    bool count(const key_type key) const {
        return occ_map_[key.int_id()];
    }
    
    template<typename... Args>
    std::pair<iterator<>, bool> emplace(key_type key, Args &&... args) {
        size_t id = key.int_id();
        if (occ_map_[id])
            return { iterator<>(id, occ_map_, data_), false };
            
        occ_map_[id] = true;
        new(data_.data() + id) mapped_type(std::forward<Args>(args)...);
        return { iterator<>(id, occ_map_, data_), true};
    }

    iterator<>      begin()         { return iterator<>(0, occ_map_, data_); }
    iterator<>      end()           { return iterator<>(iterator<>::NPOS, occ_map_, data_); }
    const_iterator  begin()  const  { return const_iterator(0, occ_map_, data_); }
    const_iterator  end()    const  { return const_iterator(const_iterator::NPOS, occ_map_, data_); }
    const_iterator  cbegin() const  { return const_iterator(0, occ_map_, data_); }
    const_iterator  cend()   const  { return const_iterator(const_iterator::NPOS, occ_map_, data_); }

  private:
    std::vector<mapped_type> data_;
    std::vector<bool> occ_map_;
};  

};

