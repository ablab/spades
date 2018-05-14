#pragma once

#include "cursor.hpp"

#include "sequence/aa.hpp"

#include <llvm/ADT/SmallVector.h>

template <class GraphCursor>
class AAGraphCursor;

template <class GraphCursor>
auto make_aa_cursors(const std::vector<GraphCursor> &cursors);

template <class GraphCursor>
inline std::ostream &operator<<(std::ostream &os, const AAGraphCursor<GraphCursor> &cursor);

template <class GraphCursor>
class AAGraphCursor {
  using This = AAGraphCursor<GraphCursor>;

 public:
  char letter() const { return to_one_letter(aa::to_aa(c0_.letter(), c1_.letter(), c2_.letter())); }

  AAGraphCursor() = default;
  AAGraphCursor(const GraphCursor &c0, const GraphCursor &c1, const GraphCursor &c2) : c0_{c0}, c1_{c1}, c2_{c2} {}
  ~AAGraphCursor() noexcept = default;
  AAGraphCursor(const This &) = default;
  AAGraphCursor(This &&) = default;
  This &operator=(This &&) = default;
  This &operator=(const This &) = default;

  bool operator==(const AAGraphCursor &other) const { return c0_ == other.c0_ && c1_ == other.c1_ && c2_ == other.c2_; }

  bool is_empty() const { return c0_.is_empty() || c1_.is_empty() || c2_.is_empty(); }

  using EdgeId = std::decay_t<decltype(GraphCursor().edge())>;

  std::vector<This> prev() const;  // TODO implement it

  std::vector<This> next() const { return from_bases(c2_.next()); }

  std::vector<GraphCursor> nucl_cursors() const { return {c0_, c1_, c2_}; }

 private:
  GraphCursor c0_, c1_, c2_;
  friend struct std::hash<This>;
  friend auto make_aa_cursors<GraphCursor>(const std::vector<GraphCursor> &cursors);
  friend std::ostream &operator<<<GraphCursor>(std::ostream &os, const AAGraphCursor<GraphCursor> &cursor);

  static std::vector<This> from_bases(const std::vector<GraphCursor> &cursors) {
    llvm::SmallVector<llvm::SmallVector<GraphCursor, 2>, 16> nexts;
    // nexts.reserve(16);

    for (const auto &cursor : cursors) {
      for (const auto &n : cursor.next()) {
        nexts.push_back({ cursor, n });
      }
    }

    std::vector<This> result;
    result.reserve(64);
    for (const auto &n : nexts) {
      assert(n.size() == 2);
      for (const auto &n2 : n.back().next()) {
        result.emplace_back(n[0], n[1], n2);
      }
    }

    return result;
  }
};

template <class GraphCursor>
inline std::ostream &operator<<(std::ostream &os, const AAGraphCursor<GraphCursor> &cursor) {
  if (cursor.is_empty()) {
    return os << "(@)";
  } else {
    return os << "(" << cursor.c0_ << ", " << cursor.c1_ << ", " << cursor.c2_ << ")";
  }
}

template <class GraphCursor>
auto make_aa_cursors(const std::vector<GraphCursor> &cursors) {
  return AAGraphCursor<GraphCursor>::from_bases(cursors);
}

template <class GraphCursor>
auto make_aa_cursors(const GraphCursor &cursor) {
  return make_aa_cursors({cursor});
}

namespace std {
template <class GraphCursor>
struct hash<AAGraphCursor<GraphCursor>> {
  std::size_t operator()(const AAGraphCursor<GraphCursor> &p) const {
    auto h = [](const GraphCursor &c) { return std::hash<GraphCursor>()(c); };
    return std::hash<size_t>()(
        hash_size_t_pair(h(p.c0_), hash_size_t_pair(h(p.c1_), h(p.c2_))));  // TODO implement a proper hash for tuples
  }
};
}  // namespace std
