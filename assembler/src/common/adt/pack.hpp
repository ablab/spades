//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <algorithm>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#include <cassert>

#include "utils/verify.hpp"

namespace adt {

class pack {
private:
    template <typename T>
    static constexpr bool is_proper_type_v = (std::is_arithmetic<T>::value || std::is_class<T>::value ||
                                              std::is_union<T>::value || std::is_enum<T>::value) &&
                                             !std::is_volatile<T>::value;

public:
    template <typename T>
    T &acquire(T *p) {
        return acquire("", p);
    }

    template <typename T>
    T &acquire(const std::string &key, T *p) {
        static_assert(is_proper_type_v<T>, "T is not a proper type for pack structure");
        return get_or_create_storage_unit_<T>().acquire(key, p);
    }

    template <typename T>
    auto &add(T &&object) {
        return add("", std::forward<T>(object));
    }

    template <typename T>
    auto &add(const std::string &key, T &&object) {
        return acquire(key, new std::remove_cv_t<std::remove_reference_t<T>>(std::forward<T>(object)));
    }

    template <typename T, typename... Args>
    T &emplace(Args &&... args) {
        return acquire("", new T(std::forward<Args>(args)...));
    }

    template <typename T, typename... Args>
    T &emplace_with_key(const std::string &key, Args &&... args) {
        return acquire(key, new T(std::forward<Args>(args)...));
    }

    template <typename T>
    size_t size() const {
        static_assert(is_proper_type_v<T>, "T is not a proper type for pack structure");
        return size(typeid(T));
    }

    size_t size() const {
        size_t size = 0;
        for (const auto &kv : storage_) {
            size += kv.second.size();
        }
        return size;
    }

    template <typename T>
    const T &get(const std::string &key = "") const {
        const auto &storage_unit = get_storage_unit_<T>();
        return storage_unit.template get_const<T>(key);
    }

    template <typename T>
    T &get_mutable(const std::string &key = "") {
        auto &storage_unit = get_storage_unit_<T>();
        return storage_unit.template get_mutable<T>(key);
    }

    template <typename T>
    T *release(const std::string &key = "") {
        auto &storage_unit = get_storage_unit_<T>();
        T *p = storage_unit.template release<T>(key);
        if (storage_unit.empty()) {
            // Remove empty storage unit
            erase(typeid(T));
        }
        return p;
    }

    bool empty() const { return storage_.empty(); }

    size_t clear() { return erase(); }

    template <typename T>
    size_t count(const std::string &key = "") const {
        static_assert(is_proper_type_v<T>, "T is not a proper type for pack structure");
        return count(typeid(T), key);
    }

    size_t erase() {
        size_t n = size();
        storage_.clear();
        return n;
    }

    template <typename T>
    size_t erase() {
        static_assert(is_proper_type_v<T>, "T is not a proper type for pack structure");
        return erase(typeid(T));
    }

    template <typename T>
    size_t erase(const std::string &key) {
        static_assert(is_proper_type_v<T>, "T is not a proper type for pack structure");
        return erase(typeid(T), key);
    }

    template <typename T>
    bool &invalidated(const std::string &key = "") {
        static_assert(is_proper_type_v<T>, "T is not a proper type for pack structure");
        return invalidated(typeid(T), key);
    }

    template <typename T>
    bool invalidated(const std::string &key = "") const {
        return const_cast<pack *>(this)->invalidated<T>(key);
    }

    template <typename T>
    void reset_invalidated(bool invalidated = false) {
        static_assert(is_proper_type_v<T>, "T is not a proper type for pack structure");
        reset_invalidated(typeid(T), invalidated);
    }

    void reset_invalidated(bool invalidated = false) {
        for (auto &kv : storage_) {
            kv.second.reset_invalidated(invalidated);
        }
    }

    size_t ntypes() const { return storage_.size(); }

    pack() = default;
    pack(pack &&) noexcept = default;
    pack &operator=(pack &&) noexcept = default;
    ~pack() noexcept = default;
    pack(const pack &) = delete;
    pack &operator=(const pack &) = delete;

protected:
    size_t size(const std::type_info &type) const {
        auto tindex = std::type_index(type);
        return storage_.count(tindex) ? storage_.find(tindex)->second.size() : 0;
    }

    size_t count(const std::type_info &type, const std::string &key = "") const {
        auto tindex = std::type_index(type);
        auto it = storage_.find(tindex);
        return it == storage_.cend() ? 0 : it->second.count(key);
    }

    size_t erase(const std::type_info &type) {
        size_t n = count(type);
        auto tindex = std::type_index(type);
        storage_.erase(tindex);
        return n;
    }

    size_t erase(const std::type_info &type, const std::string &key) {
        auto tindex = std::type_index(type);
        auto it = storage_.find(tindex);
        if (it == storage_.cend()) {
            return 0;
        }
        size_t n = it->second.erase(key);
        if (it->second.empty()) {
            // Remove empty storage unit
            erase(type);
        }
        return n;
    }

    bool &invalidated(const std::type_info &type, const std::string &key = "") {
        VERIFY(count(type, key));
        auto tindex = std::type_index(type);
        return storage_.find(tindex)->second.invalidated(key);
    }

    void reset_invalidated(const std::type_info &type, bool invalidated = false) {
        auto tindex = std::type_index(type);
        auto it = storage_.find(tindex);
        if (it != storage_.cend()) {
            it->second.reset_invalidated(invalidated);
        }
    }

private:
    class storage_unit {
    private:
        template <typename T>
        static void destroy(void *p) {
            delete static_cast<T *>(p);
        }

        std::type_index type_;
        void (*destroy_)(void *);

        struct record {
            void *pointer;
            bool invalidated;
        };

        std::unordered_map<std::string, record> records_;

        template <typename T>
        storage_unit(T *) : type_{std::type_index(typeid(T))}, destroy_{&destroy<T>} {}

    public:
        ~storage_unit() noexcept {
            for (auto &kv : records_) {
                destroy_(kv.second.pointer);
            }
        }

        storage_unit(storage_unit &&) noexcept = default;
        storage_unit &operator=(storage_unit &&) noexcept = default;
        storage_unit() = delete;
        storage_unit(const storage_unit &) = delete;
        storage_unit &operator=(const storage_unit &) = delete;

        template <typename T>
        static storage_unit make_storage_unit() {
            return storage_unit(static_cast<T *>(nullptr));
        }

        size_t size() const { return records_.size(); }

        bool empty() const { return records_.empty(); }

        size_t count(const std::string &key) const { return records_.count(key); }

        auto get_iterator_(const std::string &key) {
            auto it = records_.find(key);
            VERIFY(it != records_.end());
            return it;
        }

        auto get_iterator_(const std::string &key) const {
            auto it = records_.find(key);
            VERIFY(it != records_.cend());
            return it;
        }

        template <typename T>
        T *release(const std::string &key) {
            VERIFY(std::type_index(typeid(T)) == type_);
            auto it = get_iterator_(key);
            T *p = static_cast<T *>(it->second.pointer);
            records_.erase(it);
            return p;
        }

        template <typename T>
        const T &get_const(const std::string &key) const {
            VERIFY(std::type_index(typeid(T)) == type_);
            auto it = get_iterator_(key);
            return *static_cast<const T *>(it->second.pointer);
        }

        template <typename T>
        T &get_mutable(const std::string &key) {
            VERIFY(std::type_index(typeid(T)) == type_);
            auto it = get_iterator_(key);
            it->second.invalidated = true;
            return *static_cast<T *>(it->second.pointer);
        }

        size_t erase(const std::string &key) {
            auto it = records_.find(key);
            if (it == records_.end()) {
                return 0;
            } else {
                destroy_(it->second.pointer);
                records_.erase(it);
                return 1;
            }
        }

        template <typename T>
        T &acquire(const std::string &key, T *p) {
            VERIFY(std::type_index(typeid(T)) == type_);
            VERIFY(!records_.count(key));
            // We should not use empty and non-empty keys at the same time
            VERIFY(!records_.count(""));
            VERIFY(!key.empty() || records_.empty());
            records_[key] = {p, true};
            return *p;
        }

        bool &invalidated(const std::string &key) {
            auto it = get_iterator_(key);
            return it->second.invalidated;
        }

        void reset_invalidated(bool invalidated = false) {
            for (auto &kv : records_) {
                kv.second.invalidated = invalidated;
            }
        }
    };

    std::unordered_map<std::type_index, storage_unit> storage_;

    template <typename T>
    storage_unit &get_or_create_storage_unit_() {
        auto tindex = std::type_index(typeid(T));
        auto it_fl = storage_.insert(std::make_pair(tindex, storage_unit::make_storage_unit<T>()));
        return it_fl.first->second;
    }

    template <typename T>
    storage_unit &get_storage_unit_() {
        static_assert(is_proper_type_v<T>, "T is not a proper type for pack structure");
        auto tindex = std::type_index(typeid(T));
        auto it = storage_.find(tindex);
        VERIFY(it != storage_.end());
        return it->second;
    }

    template <typename T>
    const storage_unit &get_storage_unit_() const {
        return const_cast<pack *>(this)->get_storage_unit_<T>();
    }
};

}  // namespace adt

// vim: set ts=4 sw=4 et :
