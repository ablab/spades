#pragma once

#include "gqf/gqf.h"
#include <mutex>
#include <cmath>
#include <cstring>
#include <functional>

namespace qf {

class cqf {
    cqf(const cqf&) = delete;
    cqf& operator=(const cqf&) = delete;

  public:
    /// The hash digest type.
    typedef uint64_t digest;

    ~cqf() { qf_destroy(&qf_); }

    cqf(uint64_t maxn)
            : insertions_(0) {
        unsigned qbits = unsigned(ceil(log2(double(maxn)))) + 1;
        num_hash_bits_ = qbits + 8;
        num_slots_ = (1ULL << qbits);
        qf_init(&qf_, num_slots_, num_hash_bits_, 0, 42);
        range_mask_ = qf_.metadata->range - 1;
        assert((range_mask_ & qf_.metadata->range) == 0);
        // fprintf(stderr, "%llu %llu %u %llu\n", maxn, num_slots_, num_hash_bits_, qf_.metadata->range);
    }
    cqf(uint64_t num_slots, unsigned hash_bits)
            : num_hash_bits_(hash_bits), num_slots_(num_slots), insertions_(0) {
        qf_init(&qf_, num_slots_, num_hash_bits_, 0, 239);
        range_mask_ = qf_.metadata->range - 1;
        assert((range_mask_ & qf_.metadata->range) == 0);
        // fprintf(stderr, "%llu %u %llu\n", num_slots_, num_hash_bits_, qf_.metadata->range);
    }

    cqf(cqf&&) noexcept = default;

    bool add(digest d, uint64_t count = 1,
             bool lock = true, bool spin = true) {
        bool res = qf_insert(&qf_, d & range_mask_, 0, count, lock, spin);
        if (res)
            insertions_ += 1;
        return res;
    }

    void expand() {
        // Create new QF having the same hash bits, but double the slots
        QF nqf;
        qf_init(&nqf, 2 * num_slots_, num_hash_bits_, 0, 239);

        merge(&nqf, &qf_);

        qf_destroy(&qf_);
        memcpy(&qf_, &nqf, sizeof(qf_));
        num_slots_ = 2 * num_slots_;
    }

    void merge(cqf &other) {
        merge(&qf_, &other.qf_);
        other.clear();
    }

    void clear() {
        qf_reset(&qf_);
        insertions_ = 0;
    }

    bool full() const {
        return occupied_slots() >= uint64_t(0.95 * double(slots()));
    }

    size_t insertions() const { return insertions_; }
    unsigned hash_bits() const { return num_hash_bits_; }
    uint64_t range_mask() const { return range_mask_; }
    uint64_t slots() const { return num_slots_; }
    uint64_t occupied_slots() const { return qf_.metadata->noccupied_slots; }
    uint64_t distinct() const { return qf_.metadata->ndistinct_elts; }

    size_t lookup(digest d, bool lock = false) const {
        return qf_count_key_value(&qf_, d & range_mask_, 0, lock);
    }

private:
    void merge(QF *qf, QF *other) {
        QFi other_cfi;

        if (qf_iterator(other, &other_cfi, 0)) {
            do {
                uint64_t key = 0, value = 0, count = 0;
                qfi_get(&other_cfi, &key, &value, &count);
                qf_insert(qf, key, value, count, true, true);
            } while (!qfi_next(&other_cfi));
        }
    }

    QF qf_;
    unsigned num_hash_bits_;
    uint64_t num_slots_;
    size_t insertions_;
    uint64_t range_mask_;
};

template<class T>
class cqf_with_hasher : public cqf {
  public:
    using cqf::digest;
    using cqf::add;
    using cqf::lookup;

    /// The hash function type.
    typedef std::function<digest(const T&)> hasher;

    cqf_with_hasher(uint64_t maxn, hasher h = nullptr)
            : cqf(maxn), hasher_(std::move(h)) {}

    cqf_with_hasher(uint64_t num_slots, unsigned hash_bits, hasher h = nullptr)
            : hasher_(std::move(h)), cqf(num_slots, hash_bits) {}

    void replace_hasher(hasher h) {
        hasher_ = std::move(h);
        clear();
    }

    bool add(const T &o, uint64_t count = 1,
             bool lock = true, bool spin = true) {
        digest d = hasher_(o);
        return add(d, count, lock, spin);
    }

    size_t lookup(const T &o, bool lock = false) const {
        digest d = hasher_(o);
        return lookup(d, lock);
    }

private:
    hasher hasher_;
};


};
