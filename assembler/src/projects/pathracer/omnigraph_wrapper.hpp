#pragma once

#include "cursor.hpp"
#include "pathtree.hpp"
#include "utils.hpp"
#include "aa_cursor.hpp"
#include "string_cursor.hpp"
#include "cached_cursor.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/components/graph_component.hpp"

#include <vector>

struct IdHolder {
    uint64_t id_ : 26;
    uint64_t pos_ : 28;

    IdHolder()
            : id_(0), pos_(-1) {}

    IdHolder(debruijn_graph::EdgeId id, uint64_t pos)
            : id_(id.int_id()), pos_(pos) {}

    explicit operator bool() const { return id_; }

    bool operator==(const IdHolder &other) const {
        return (id_ == other.id_) && (pos_ == other.pos_);
    }

    debruijn_graph::EdgeId e() const {
        return id_;
    }

    size_t pos() const {
        return pos_;
    };
};

static_assert(sizeof(IdHolder) == sizeof(uint64_t), "Invalid packed id holder");

class DebruijnGraphCursor {
  public:
      DebruijnGraphCursor()
              : holder_{} {}

      DebruijnGraphCursor(debruijn_graph::EdgeId e, size_t position) : holder_{e, position} {}

      DebruijnGraphCursor(const DebruijnGraphCursor &) = default;
      DebruijnGraphCursor(DebruijnGraphCursor &&) = default;
      DebruijnGraphCursor &operator=(const DebruijnGraphCursor &) = default;
      DebruijnGraphCursor &operator=(DebruijnGraphCursor &&) = default;
      ~DebruijnGraphCursor() noexcept = default;

      bool operator==(const DebruijnGraphCursor &other) const {
          return holder_ == other.holder_;
    }

    using EdgeId = debruijn_graph::EdgeId;
    EdgeId edge() const { return holder_.e(); }
    size_t position() const { return holder_.pos(); }
    bool is_empty() const { return !holder_; }

    char letter(const void *context) const {
        return nucl(g(context).EdgeNucls(edge())[position()]); }

    // bool is_convergent(const void *context) const {
    //     const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);
    //
    //     return position() == g.k() &&
    //             g.IncomingEdgeCount(g.EdgeStart(edge())) > 1;
    //     // return prev().size() > 1;
    // }
    //
    // bool is_divergent(const void *context) const {
    //     const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);
    //
    //     return position() + 1 == g.length(edge()) + g.k() &&
    //             g.OutgoingEdgeCount(g.EdgeEnd(edge())) > 1;
    //     // return next().size() > 1;
    // }

    std::vector<DebruijnGraphCursor> prev(const void*) const;
    std::vector<DebruijnGraphCursor> next(const void*) const;

    static std::vector<DebruijnGraphCursor> get_cursors(const debruijn_graph::ConjugateDeBruijnGraph &g, const debruijn_graph::EdgeId &e, size_t pos) {
        // Unfortunately, several different cursors (actually different, having different prev's) may have the same edge & position
        // Therefore, it's impossible to design a correct get_cursor() method, we have to implement get_cursorS()
        DebruijnGraphCursor cursor{e, pos};
        std::vector<DebruijnGraphCursor> result;
        cursor.generate_normalized_cursors(result, &g);
        return result;
    }

    static std::vector<DebruijnGraphCursor> all(const debruijn_graph::ConjugateDeBruijnGraph &g);

  private:
    const debruijn_graph::ConjugateDeBruijnGraph &g(const void *context) const {
        return *static_cast<const debruijn_graph::ConjugateDeBruijnGraph *>(context);
    }
    
    void normalize_prefix_to_suffix(const void *context) {
        const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);

        // This method is used ONLY in all() generators for duplicates merging
        // normalization is not one-valued, but for duplicates deletion it's sufficient to normalize to any of normal variants
        // Also please not that in terms of next() all normalized versions are equivalent
        // Thus, for initial cursors we can (and should!) perform even more aggressive normalization "to brother":
        // if we are in vertex, replace the edge with the edge having the least id among all ingoing edges.
        while (position() < g.k() && g.IncomingEdgeCount(g.EdgeStart(edge())) > 0) {
            holder_ = { *g.in_begin(g.EdgeStart(edge())), g.length(edge()) + position() };
        }
    }

    void generate_normalized_cursors(std::vector<DebruijnGraphCursor> &out,
                                     const void *context) const {
        const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);
        if (position() < g.k() && g.IncomingEdgeCount(g.EdgeStart(edge())) > 0) {
            for (EdgeId new_e : g.IncomingEdges(g.EdgeStart(edge()))) {
                size_t new_position = g.length(new_e) + position();
                DebruijnGraphCursor cursor{new_e, new_position};
                cursor.generate_normalized_cursors(out, context);
            }
        } else {
            out.push_back(*this);
        }
    }

    friend struct std::hash<DebruijnGraphCursor>;
    friend std::ostream &operator<<(std::ostream &os, const DebruijnGraphCursor &p);
    
    IdHolder holder_;
};

static_assert(sizeof(DebruijnGraphCursor) == sizeof(uint64_t), "Invalid packed id holder");


class DebruijnComponentCursor {
  public:
    DebruijnComponentCursor()
            : holder_() {}

    DebruijnComponentCursor(debruijn_graph::EdgeId e, size_t position)
            : holder_{e, position} {}

    DebruijnComponentCursor(const DebruijnComponentCursor &) = default;
    DebruijnComponentCursor(DebruijnComponentCursor &&) = default;
    DebruijnComponentCursor &operator=(const DebruijnComponentCursor &) = default;
    DebruijnComponentCursor &operator=(DebruijnComponentCursor &&) = default;
    ~DebruijnComponentCursor() noexcept = default;

    bool operator==(const DebruijnComponentCursor &other) const {
        return holder_ == other.holder_;
    }

    using EdgeId = debruijn_graph::EdgeId;
    EdgeId edge() const { return holder_.e(); }
    size_t position() const { return holder_.pos(); }

    bool is_empty() const { return !holder_; }

    char letter(const void *context) const { return nucl(this->g(context).EdgeNucls(edge())[position()]); }

    // bool is_convergent(const void *context) const {
    //     const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);
    //     return position() == g.k() &&
    //             (g.IncomingEdgeCount(g.EdgeStart(edge())) > 1 && !this->c(context).IsBorder(g.EdgeStart(edge())));
    //     // return prev().size() > 1;
    // }
    //
    // bool is_divergent(const void *context) const {
    //     const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);
    //     return position() + 1 == g.length(edge()) + g.k() &&
    //             (g.OutgoingEdgeCount(g.EdgeEnd(edge())) > 1 && !this->c(context).IsBorder(g.EdgeEnd(edge())));
    //     // return next().size() > 1;
    // }

    std::vector<DebruijnComponentCursor> prev(const void*) const;
    std::vector<DebruijnComponentCursor> next(const void*) const;

    static DebruijnComponentCursor get_cursor(const omnigraph::GraphComponent<debruijn_graph::ConjugateDeBruijnGraph> &c,
                                              const debruijn_graph::EdgeId &e, size_t pos) {
        DebruijnComponentCursor cursor{e, pos};
        cursor.normalize_prefix_to_suffix(&c);
        return cursor;
    }

    static std::vector<DebruijnComponentCursor> all(const omnigraph::GraphComponent<debruijn_graph::ConjugateDeBruijnGraph> &c);

  private:
    const omnigraph::GraphComponent<debruijn_graph::ConjugateDeBruijnGraph> &c(const void *context) const {
        return *static_cast<const omnigraph::GraphComponent<debruijn_graph::ConjugateDeBruijnGraph> *>(context);
    }
    
    const debruijn_graph::ConjugateDeBruijnGraph &g(const void *context) const {
        return this->c(context).g();
    }

    void normalize_prefix_to_suffix(const void *context) {
        const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);
        while (position() < g.k() &&
               (g.IncomingEdgeCount(g.EdgeStart(edge())) > 0 && !this->c(context).IsBorder(g.EdgeStart(edge())))) {
            holder_ = { *g.in_begin(g.EdgeStart(edge())), g.length(edge()) + position() };
        }
    }

    friend struct std::hash<DebruijnComponentCursor>;
    friend std::ostream &operator<<(std::ostream &os, const DebruijnComponentCursor &p);
        
    IdHolder holder_;
};

static_assert(sizeof(DebruijnComponentCursor) == sizeof(uint64_t), "Invalid packed id holder");


namespace std {
template <>
struct hash<DebruijnGraphCursor> {
    std::size_t operator()(const DebruijnGraphCursor &p) const {
        return std::hash<size_t>()(hash_size_t_pair(p.edge().hash(), p.position()));
    }
};
}

namespace std {
template <>
struct hash<DebruijnComponentCursor> {
    std::size_t operator()(const DebruijnComponentCursor &p) const {
        return std::hash<size_t>()(hash_size_t_pair(p.edge().hash(), p.position()));
    }
};
}

inline std::ostream &operator<<(std::ostream &os, const DebruijnGraphCursor &p) {
  if (p.is_empty()) {
    return os << "(@)";
  } else {
    return os << "(" << p.edge() << ", " << p.position() << ")";
  }
}

inline std::ostream &operator<<(std::ostream &os, const DebruijnComponentCursor &p) {
  if (p.is_empty()) {
    return os << "(@)";
  } else {
    return os << "(" << p.edge() << ", " << p.position() << ")";
  }
}

namespace hmm {
struct Fees;
};

PathSet<DebruijnGraphCursor> find_best_path(const hmm::Fees &fees,
                                            const std::vector<DebruijnGraphCursor> &initial,
                                            const void *context);
PathSet<DebruijnComponentCursor> find_best_path(const hmm::Fees &fees,
                                                const std::vector<DebruijnComponentCursor> &initial,
                                                const void *context);
PathSet<ReversalGraphCursor<DebruijnGraphCursor>> find_best_path_rev(const hmm::Fees &fees,
                                                                     const std::vector<ReversalGraphCursor<DebruijnGraphCursor>> &initial,
                                                                     const void *context);
PathSet<ReversalGraphCursor<DebruijnComponentCursor>> find_best_path_rev(const hmm::Fees &fees,
                                                                         const std::vector<ReversalGraphCursor<DebruijnComponentCursor>> &initial,
                                                                         const void *context);
PathSet<AAGraphCursor<DebruijnComponentCursor>> find_best_path(const hmm::Fees &fees,
                                                               const std::vector<AAGraphCursor<DebruijnComponentCursor>> &initial,
                                                               const void *context);
PathSet<AAGraphCursor<DebruijnGraphCursor>> find_best_path(const hmm::Fees &fees,
                                                           const std::vector<AAGraphCursor<DebruijnGraphCursor>> &initial,
                                                           const void *context);
PathSet<AAGraphCursor<RestrictedGraphCursor<DebruijnGraphCursor>>> find_best_path(const hmm::Fees &fees,
                                                                                  const std::vector<AAGraphCursor<RestrictedGraphCursor<DebruijnGraphCursor>>> &initial,
                                                                                  const void *context);
PathSet<RestrictedGraphCursor<DebruijnGraphCursor>> find_best_path(const hmm::Fees &fees,
                                                                   const std::vector<RestrictedGraphCursor<DebruijnGraphCursor>> &initial,
                                                                   const void *context);
PathSet<AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>> find_best_path(const hmm::Fees &fees,
                                                                                           const std::vector<AAGraphCursor<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>>> &initial,
                                                                                           const void *context);
PathSet<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>> find_best_path(const hmm::Fees &fees,
                                                                            const std::vector<OptimizedRestrictedGraphCursor<DebruijnGraphCursor>> &initial,
                                                                            const void *context);
PathSet<StringCursor> find_best_path(const hmm::Fees &fees,
                                     const std::vector<StringCursor> &initial,
                                     const void *context);
PathSet<CachedCursor> find_best_path(const hmm::Fees &fees,
                                     const std::vector<CachedCursor> &initial,
                                     const void *context);
