#pragma once

#include "gqf/gqf.h"

namespace qf {

template<class T>
class cqf {
    cqf(const cqf&) = delete;
    cqf& operator=(const cqf&) = delete;

  public:
    /// The hash digest type.
    typedef uint64_t digest;

    /// The hash function type.
    typedef std::function<digest(const T&, uint64_t seed)> hasher;

    ~cqf() { qf_destroy(&qf_); }

    cqf(hasher h, uint64_t maxn)
            : hasher_(std::move(h)), insertions_(0) {
        unsigned qbits = unsigned(ceil(log2(maxn))) + 1;
        num_hash_bits_ = qbits + 8;
        num_slots_ = (1ULL << qbits);
        qf_init(&qf_, num_slots_, num_hash_bits_, 0, 42);
        fprintf(stderr, "%llu %llu %u %llu\n", maxn, num_slots_, num_hash_bits_, qf_.metadata->range);
    }
    cqf(hasher h, uint64_t num_slots, unsigned hash_bits)
            : hasher_(std::move(h)),
              num_hash_bits_(hash_bits), num_slots_(num_slots), insertions_(0) {
        qf_init(&qf_, num_slots_, num_hash_bits_, 0, 239);
        fprintf(stderr, "%llu %u %llu\n", num_slots_, num_hash_bits_, qf_.metadata->range);
    }

    cqf(cqf&&) = default;

    bool add(digest d, uint64_t count = 1,
             bool lock = true, bool spin = true) {
        return qf_insert(&qf_, d % qf_.metadata->range, 0, count, lock, spin);
    }

    bool add(const T &o, uint64_t count = 1,
             bool lock = true, bool spin = true) {
        digest d = hasher_(o, 0);
        return add(d, count, lock, spin);
    }

    void merge(cqf<T> &other) {
        QFi other_cfi;

        if (qf_iterator(&other.qf_, &other_cfi, 0)) {
            do {
                uint64_t key = 0, value = 0, count = 0;
                qfi_get(&other_cfi, &key, &value, &count);
                qf_insert(&qf_, key, value, count, true, true);
            } while (!qfi_next(&other_cfi));
        }
        other.clear();
    }

    void clear() {
        qf_reset(&qf_);
        insertions_ = 0;
    }

    size_t insertions() const { return insertions_; }

    unsigned hash_bits() const {
        return num_hash_bits_;
    }

    uint64_t slots() const {
        return num_slots_;
    }

    size_t lookup(digest d) const {
        return qf_count_key_value(&qf_, d % qf_.metadata->range, 0);
    }

    size_t lookup(const T&o) const {
        digest d = hasher_(o, 0);
        return lookup(d);
    }

protected:
    hasher hasher_;
    QF qf_;
    unsigned num_hash_bits_;
    uint64_t num_slots_;
    size_t insertions_;
};


};
