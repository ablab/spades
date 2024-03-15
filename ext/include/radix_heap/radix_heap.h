#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace radix_heap {
namespace internal {
template<bool Is64bit> class find_bucket_impl;

template<>
class find_bucket_impl<false> {
 public:
  static inline constexpr size_t find_bucket(uint32_t x, uint32_t last) {
    return x == last ? 0 : 32 - __builtin_clz(x ^ last);
  }
};

template<>
class find_bucket_impl<true> {
 public:
  static inline constexpr size_t find_bucket(uint64_t x, uint64_t last) {
    return x == last ? 0 : 64 - __builtin_clzll(x ^ last);
  }
};

template<typename T>
inline constexpr size_t find_bucket(T x, T last) {
  return find_bucket_impl<sizeof(T) == 8>::find_bucket(x, last);
}

template<typename KeyType, bool IsSigned> class encoder_impl_integer;

template<typename KeyType>
class encoder_impl_integer<KeyType, false> {
 public:
  typedef KeyType key_type;
  typedef KeyType unsigned_key_type;

  inline static constexpr unsigned_key_type encode(key_type x) {
    return x;
  }

  inline static constexpr key_type decode(unsigned_key_type x) {
    return x;
  }
};

template<typename KeyType>
class encoder_impl_integer<KeyType, true> {
 public:
  typedef KeyType key_type;
  typedef typename std::make_unsigned<KeyType>::type unsigned_key_type;

  inline static constexpr unsigned_key_type encode(key_type x) {
    return static_cast<unsigned_key_type>(x) ^
        (unsigned_key_type(1) << unsigned_key_type(std::numeric_limits<unsigned_key_type>::digits - 1));
  }

  inline static constexpr key_type decode(unsigned_key_type x) {
    return static_cast<key_type>
        (x ^ (unsigned_key_type(1) << (std::numeric_limits<unsigned_key_type>::digits - 1)));
  }
};

template<typename KeyType, typename UnsignedKeyType>
class encoder_impl_decimal {
public:
  typedef KeyType key_type;
  typedef UnsignedKeyType unsigned_key_type;

  inline static constexpr unsigned_key_type encode(key_type x) {
    return raw_cast<key_type, unsigned_key_type>(x) ^
        ((-(raw_cast<key_type, unsigned_key_type>(x) >> (std::numeric_limits<unsigned_key_type>::digits - 1))) |
         (unsigned_key_type(1) << (std::numeric_limits<unsigned_key_type>::digits - 1)));
  }

  inline static constexpr key_type decode(unsigned_key_type x) {
    return raw_cast<unsigned_key_type, key_type>
        (x ^ (((x >> (std::numeric_limits<unsigned_key_type>::digits - 1)) - 1) |
              (unsigned_key_type(1) << (std::numeric_limits<unsigned_key_type>::digits - 1))));
  }

 private:
  template<typename T, typename U>
  union raw_cast {
   public:
    constexpr raw_cast(T t) : t_(t) {}
    operator U() const { return u_; }

   private:
    T t_;
    U u_;
  };
};

template<typename KeyType>
class encoder : public encoder_impl_integer<KeyType, std::is_signed<KeyType>::value> {};
template<>
class encoder<float> : public encoder_impl_decimal<float, uint32_t> {};
template<>
class encoder<double> : public encoder_impl_decimal<double, uint64_t> {};
}  // namespace internal

template<typename KeyType, typename EncoderType = internal::encoder<KeyType>>
class radix_heap {
 public:
  typedef KeyType key_type;
  typedef EncoderType encoder_type;
  typedef typename encoder_type::unsigned_key_type unsigned_key_type;

  radix_heap() : size_(0), last_(), buckets_() {
    buckets_min_.fill(std::numeric_limits<unsigned_key_type>::max());
  }

  void push(key_type key) {
    const unsigned_key_type x = encoder_type::encode(key);
    assert(last_ <= x);
    ++size_;
    const size_t k = internal::find_bucket(x, last_);
    buckets_[k].emplace_back(x);
    buckets_min_[k] = std::min(buckets_min_[k], x);
  }

  key_type top() {
    pull();
    return encoder_type::decode(last_);
  }

  void pop() {
    pull();
    buckets_[0].pop_back();
    --size_;
  }

  size_t size() const {
    return size_;
  }

  bool empty() const {
    return size_ == 0;
  }

  void clear() {
    size_ = 0;
    last_ = key_type();
    for (auto &b : buckets_) b.clear();
    buckets_min_.fill(std::numeric_limits<unsigned_key_type>::max());
  }

  void swap(radix_heap<KeyType, EncoderType> &a) {
    std::swap(size_, a.size_);
    std::swap(last_, a.last_);
    buckets_.swap(a.buckets_);
    buckets_min_.swap(a.buckets_min_);
  }

 private:
  size_t size_;
  unsigned_key_type last_;
  std::array<std::vector<unsigned_key_type>,
             std::numeric_limits<unsigned_key_type>::digits + 1> buckets_;
  std::array<unsigned_key_type,
             std::numeric_limits<unsigned_key_type>::digits + 1> buckets_min_;

  void pull() {
    assert(size_ > 0);
    if (!buckets_[0].empty()) return;

    size_t i;
    for (i = 1; buckets_[i].empty(); ++i);
    last_ = buckets_min_[i];

    for (unsigned_key_type x : buckets_[i]) {
      const size_t k = internal::find_bucket(x, last_);
      buckets_[k].emplace_back(x);
      buckets_min_[k] = std::min(buckets_min_[k], x);
    }
    buckets_[i].clear();
    buckets_min_[i] = std::numeric_limits<unsigned_key_type>::max();
  }
};

template<typename KeyType, typename ValueType, typename EncoderType = internal::encoder<KeyType>>
class pair_radix_heap {
 public:
  typedef KeyType key_type;
  typedef ValueType value_type;
  typedef EncoderType encoder_type;
  typedef typename encoder_type::unsigned_key_type unsigned_key_type;

  pair_radix_heap() : size_(0), last_(), buckets_() {
    buckets_min_.fill(std::numeric_limits<unsigned_key_type>::max());
  }

  void push(key_type key, const value_type &value) {
    const unsigned_key_type x = encoder_type::encode(key);
    assert(last_ <= x);
    ++size_;
    const size_t k = internal::find_bucket(x, last_);
    buckets_[k].emplace_back(x, value);
    buckets_min_[k] = std::min(buckets_min_[k], x);
  }

  void push(key_type key, value_type &&value) {
    const unsigned_key_type x = encoder_type::encode(key);
    assert(last_ <= x);
    ++size_;
    const size_t k = internal::find_bucket(x, last_);
    buckets_[k].emplace_back(x, std::move(value));
    buckets_min_[k] = std::min(buckets_min_[k], x);
  }

  template <class... Args>
  void emplace(key_type key, Args&&... args) {
    const unsigned_key_type x = encoder_type::encode(key);
    assert(last_ <= x);
    ++size_;
    const size_t k = internal::find_bucket(x, last_);
    buckets_[k].emplace_back(std::piecewise_construct,
                             std::forward_as_tuple(x), std::forward_as_tuple(args...));
    buckets_min_[k] = std::min(buckets_min_[k], x);
  }

  key_type top_key() {
    pull();
    return encoder_type::decode(last_);
  }

  value_type &top_value() {
    pull();
    return buckets_[0].back().second;
  }

  void pop() {
    pull();
    buckets_[0].pop_back();
    --size_;
  }

  size_t size() const {
    return size_;
  }

  bool empty() const {
    return size_ == 0;
  }

  void clear() {
    size_ = 0;
    last_ = key_type();
    for (auto &b : buckets_) b.clear();
    buckets_min_.fill(std::numeric_limits<unsigned_key_type>::max());
  }

  void swap(pair_radix_heap<KeyType, ValueType, EncoderType> &a) {
    std::swap(size_, a.size_);
    std::swap(last_, a.last_);
    buckets_.swap(a.buckets_);
    buckets_min_.swap(a.buckets_min_);
  }

 private:
  size_t size_;
  unsigned_key_type last_;
  std::array<std::vector<std::pair<unsigned_key_type, value_type>>,
             std::numeric_limits<unsigned_key_type>::digits + 1> buckets_;
  std::array<unsigned_key_type,
             std::numeric_limits<unsigned_key_type>::digits + 1> buckets_min_;

  void pull() {
    assert(size_ > 0);
    if (!buckets_[0].empty()) return;

    size_t i;
    for (i = 1; buckets_[i].empty(); ++i);
    last_ = buckets_min_[i];

    for (size_t j = 0; j < buckets_[i].size(); ++j) {
      const unsigned_key_type x = buckets_[i][j].first;
      const size_t k = internal::find_bucket(x, last_);
      buckets_[k].emplace_back(std::move(buckets_[i][j]));
      buckets_min_[k] = std::min(buckets_min_[k], x);
    }
    buckets_[i].clear();
    buckets_min_[i] = std::numeric_limits<unsigned_key_type>::max();
  }
};
}  // namespace radix_heap
