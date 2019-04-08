#pragma once

#include "utils.hpp"

#include "assembly_graph/core/graph.hpp"

#include <cereal/archives/binary.hpp>
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

    template <class Archive>
    void serialize(Archive &archive) {
        archive(binary_pod(this));
    }
};

static_assert(sizeof(IdHolder) == sizeof(uint64_t), "Invalid packed id holder");

class DebruijnGraphCursor {
public:
    using Context = const debruijn_graph::ConjugateDeBruijnGraph*;
    DebruijnGraphCursor() : holder_{} {}

    DebruijnGraphCursor(debruijn_graph::EdgeId e, size_t position) : holder_{e, position} {}

    DebruijnGraphCursor(const DebruijnGraphCursor &) = default;
    DebruijnGraphCursor(DebruijnGraphCursor &&) = default;
    DebruijnGraphCursor &operator=(const DebruijnGraphCursor &) = default;
    DebruijnGraphCursor &operator=(DebruijnGraphCursor &&) = default;
    ~DebruijnGraphCursor() noexcept = default;

    bool operator==(const DebruijnGraphCursor &other) const { return holder_ == other.holder_; }

    using EdgeId = debruijn_graph::EdgeId;
    EdgeId edge() const { return holder_.e(); }
    size_t position() const { return holder_.pos(); }
    bool is_empty() const { return !holder_; }  //FIXME check it!!!!!

    char letter(Context context) const {
        return nucl2(g(context).EdgeNucls(edge())[position()]); }

    template <class Archive>
    void serialize(Archive &archive) {
        archive(holder_);
    }
    // bool is_convergent(Context context) const {
    //     const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);
    //
    //     return position() == g.k() &&
    //             g.IncomingEdgeCount(g.EdgeStart(edge())) > 1;
    //     // return prev().size() > 1;
    // }
    //
    // bool is_divergent(Context context) const {
    //     const debruijn_graph::ConjugateDeBruijnGraph &g = this->g(context);
    //
    //     return position() + 1 == g.length(edge()) + g.k() &&
    //             g.OutgoingEdgeCount(g.EdgeEnd(edge())) > 1;
    //     // return next().size() > 1;
    // }

    std::vector<DebruijnGraphCursor> prev(Context) const;
    std::vector<DebruijnGraphCursor> next(Context) const;

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
    const debruijn_graph::ConjugateDeBruijnGraph &g(Context context) const {
        return *context;
    }

    void normalize_prefix_to_suffix(Context context) {
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
                                     Context context) const {
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

namespace std {
template <>
struct hash<DebruijnGraphCursor> {
    std::size_t operator()(const DebruijnGraphCursor &p) const {
        return std::hash<size_t>()(hash_size_t_pair(p.edge().hash(), p.position()));
    }
};
}  // namespace std

inline std::ostream &operator<<(std::ostream &os, const DebruijnGraphCursor &p) {
  if (p.is_empty()) {
    return os << "(@)";
  } else {
    return os << "(" << p.edge() << ", " << p.position() << ")";
  }
}
inline bool operator<(const DebruijnGraphCursor &c1, const DebruijnGraphCursor &c2) {
    return std::make_tuple(c1.edge(), c1.position()) < std::make_tuple(c2.edge(), c2.position());
}
