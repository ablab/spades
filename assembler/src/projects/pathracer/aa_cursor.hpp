#pragma once

#include "utils.hpp"

#include "sequence/aa.hpp"

#include <llvm/ADT/SmallVector.h>

template <class GraphCursor>
class AAGraphCursor;

template <class GraphCursor>
auto make_aa_cursors(const std::vector<GraphCursor> &cursors, typename GraphCursor::Context context);

template <class GraphCursor>
inline std::ostream &operator<<(std::ostream &os, const AAGraphCursor<GraphCursor> &cursor);

template <class GraphCursor>
class AAGraphCursor {
  using This = AAGraphCursor<GraphCursor>;

 public:
  using Context = typename GraphCursor::Context;
  char letter(Context context) const {
      return c0_.is_empty() || c1_.is_empty() || c2_.is_empty() ?
          '*' :  // FIXME use special sign for frame shift
          to_one_letter(aa::to_aa(c0_.letter(context), c1_.letter(context), c2_.letter(context)));
  }

  template <class Archive>
  void serialize(Archive &archive) {
    archive(c0_, c1_, c2_);
  }
  AAGraphCursor() = default;
  AAGraphCursor(const GraphCursor &c0, const GraphCursor &c1, const GraphCursor &c2) : c0_{c0}, c1_{c1}, c2_{c2} {}
  ~AAGraphCursor() noexcept = default;
  AAGraphCursor(const This &) = default;
  AAGraphCursor(This &&) = default;
  This &operator=(This &&) = default;
  This &operator=(const This &) = default;

  bool operator==(const AAGraphCursor &other) const { return c0_ == other.c0_ && c1_ == other.c1_ && c2_ == other.c2_; }

  bool is_empty() const { return c0_.is_empty() && c1_.is_empty() && c2_.is_empty(); }

  using EdgeId = std::decay_t<decltype(GraphCursor().edge())>;

  std::vector<This> prev(Context context) const { return from_bases_prev(c0_.prev(context), context); }

  std::vector<This> next(Context context) const { return from_bases_next(c2_.next(context), context); }
  std::vector<This> next_frame_shift(Context context) const { return from_bases_next12(c2_.next(context), context); }

  std::vector<GraphCursor> nucl_cursors() const { return {c0_, c1_, c2_}; }

 private:
  GraphCursor c0_, c1_, c2_;
  friend struct std::hash<This>;
  friend auto make_aa_cursors<GraphCursor>(const std::vector<GraphCursor> &cursors, Context context);
  friend std::ostream &operator<<<GraphCursor>(std::ostream &os, const AAGraphCursor<GraphCursor> &cursor);

  static std::vector<This> from_bases_next12(const std::vector<GraphCursor> &cursors,
                                             Context context) {
    std::vector<This> result;
    result.reserve(20);

    for (const auto &cursor : cursors) {
      result.emplace_back(GraphCursor(), GraphCursor(), cursor);
      for (const auto &n : cursor.next(context)) {
        result.emplace_back(GraphCursor(), cursor, n);
      }
    }

    return result;
  }

  static std::vector<This> from_bases_next(const std::vector<GraphCursor> &cursors,
                                           Context context) {
    std::vector<This> result;
    result.reserve(64);

    for (const auto &cursor : cursors)
      for (const auto &n1 : cursor.next(context))
        for (const auto &n2 : n1.next(context))
          nexts.push_back({ cursor, n1, n2 });

    return result;
  }

  static std::vector<This> from_bases_prev(const std::vector<GraphCursor> &cursors,
                                           Context context) {
    std::vector<This> result;
    result.reserve(64);

    for (const auto &cursor : cursors)
      for (const auto &p1 : cursor.prev(context))
        for (const auto &p2 : p1.prev(context))
          result.emplace_back(p2, p1, cursor);

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
auto make_aa_cursors(const std::vector<GraphCursor> &cursors, typename GraphCursor::Context context) {
  return AAGraphCursor<GraphCursor>::from_bases_next(cursors, context);
}

template <class GraphCursor>
auto make_aa_cursors(const GraphCursor &cursor, typename GraphCursor::Context context) {
  return make_aa_cursors(std::vector<GraphCursor>({cursor}), context);
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

// vim: set ts=2 sw=2 et :
