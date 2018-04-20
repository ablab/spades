#pragma once

#include "cursor.hpp"
#include "pathtree.hpp"
#include "utils.hpp"
#include "aa_cursor.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/components/graph_component.hpp"

#include <vector>

class DebruijnGraphCursor : public AbstractGraphCursor<DebruijnGraphCursor> {
  public:
    DebruijnGraphCursor()
            : g_{nullptr}, e_(), position_(-1) {}

    DebruijnGraphCursor(const debruijn_graph::ConjugateDeBruijnGraph *graph,
                        debruijn_graph::EdgeId e, size_t position)
            : g_(graph), e_(e), position_(position) {
        normalize_prefix_to_suffix();
    }

    DebruijnGraphCursor(const DebruijnGraphCursor &) = default;
    DebruijnGraphCursor(DebruijnGraphCursor &&) = default;
    DebruijnGraphCursor &operator=(const DebruijnGraphCursor &) = default;
    DebruijnGraphCursor &operator=(DebruijnGraphCursor &&) = default;
    ~DebruijnGraphCursor() noexcept = default;

    bool operator==(const DebruijnGraphCursor &other) const {
        return (e_ == other.e_) && (position_ == other.position_);
    }

    using EdgeId = debruijn_graph::EdgeId;
    EdgeId edge() const { return e_; }
    std::vector<EdgeId> edges() const { return { e_ }; }

    size_t position() const { return position_; }

    bool is_empty() const { return e_.get() == nullptr; }

    char letter() const { return nucl(g_->EdgeNucls(e_)[position_]); }

    bool is_convergent() const {
        return position_ == g_->k() &&
                g_->IncomingEdgeCount(g_->EdgeStart(e_)) > 1;
        // return prev().size() > 1;
    }

    bool is_divergent() const {
        return position_ + 1 == g_->length(e_) + g_->k() &&
                g_->OutgoingEdgeCount(g_->EdgeEnd(e_)) > 1;
        // return next().size() > 1;
    }

    std::vector<DebruijnGraphCursor> prev() const;
    std::vector<DebruijnGraphCursor> next() const;

  private:
    void normalize_prefix_to_suffix() {
        while (position_ < g_->k() && g_->IncomingEdgeCount(g_->EdgeStart(e_)) > 0) {
            e_ = *g_->in_begin(g_->EdgeStart(e_));
            position_ = g_->length(e_) + position_;
        }
    }

    friend struct std::hash<DebruijnGraphCursor>;
    friend std::ostream &operator<<(std::ostream &os, const DebruijnGraphCursor &p);

    const debruijn_graph::ConjugateDeBruijnGraph *g_;
    debruijn_graph::EdgeId e_;
    size_t position_;
};

class DebruijnComponentCursor : public AbstractGraphCursor<DebruijnComponentCursor> {
  public:
    DebruijnComponentCursor()
            : c_{nullptr}, e_(), position_(-1) {}

    DebruijnComponentCursor(const omnigraph::GraphComponent<debruijn_graph::ConjugateDeBruijnGraph> *component,
                            debruijn_graph::EdgeId e, size_t position)
            : c_(component), e_(e), position_(position) {
        normalize_prefix_to_suffix();
    }

    DebruijnComponentCursor(const DebruijnComponentCursor &) = default;
    DebruijnComponentCursor(DebruijnComponentCursor &&) = default;
    DebruijnComponentCursor &operator=(const DebruijnComponentCursor &) = default;
    DebruijnComponentCursor &operator=(DebruijnComponentCursor &&) = default;
    ~DebruijnComponentCursor() noexcept = default;

    bool operator==(const DebruijnComponentCursor &other) const {
        return (e_ == other.e_) && (position_ == other.position_);
    }

    using EdgeId = debruijn_graph::EdgeId;
    EdgeId edge() const { return e_; }
    std::vector<EdgeId> edges() const { return { e_ }; }
    size_t position() const { return position_; }

    bool is_empty() const { return e_.get() == nullptr; }

    char letter() const { return nucl(c_->g().EdgeNucls(e_)[position_]); }

    bool is_convergent() const {
        const debruijn_graph::ConjugateDeBruijnGraph &g = c_->g();
        return position_ == g.k() &&
                (g.IncomingEdgeCount(g.EdgeStart(e_)) > 1 && !c_->IsBorder(g.EdgeStart(e_)));
        // return prev().size() > 1;
    }

    bool is_divergent() const {
        const debruijn_graph::ConjugateDeBruijnGraph &g = c_->g();
        return position_ + 1 == g.length(e_) + g.k() &&
                (g.OutgoingEdgeCount(g.EdgeEnd(e_)) > 1 && !c_->IsBorder(g.EdgeEnd(e_)));
        // return next().size() > 1;
    }

    std::vector<DebruijnComponentCursor> prev() const;
    std::vector<DebruijnComponentCursor> next() const;

  private:
    void normalize_prefix_to_suffix() {
        const debruijn_graph::ConjugateDeBruijnGraph &g = c_->g();
        while (position_ < g.k() &&
               (g.IncomingEdgeCount(g.EdgeStart(e_)) > 0 && !c_->IsBorder(g.EdgeStart(e_)))) {
            e_ = *g.in_begin(g.EdgeStart(e_));
            position_ = g.length(e_) + position_;
        }
    }

    friend struct std::hash<DebruijnComponentCursor>;
    friend std::ostream &operator<<(std::ostream &os, const DebruijnComponentCursor &p);

    const omnigraph::GraphComponent<debruijn_graph::ConjugateDeBruijnGraph> *c_;
    debruijn_graph::EdgeId e_;
    size_t position_;
};

namespace std {
template <>
struct hash<DebruijnGraphCursor> {
    std::size_t operator()(const DebruijnGraphCursor &p) const {
        return std::hash<size_t>()(hash_size_t_pair(p.e_.hash(), p.position_));
    }
};
}

namespace std {
template <>
struct hash<DebruijnComponentCursor> {
    std::size_t operator()(const DebruijnComponentCursor &p) const {
        return std::hash<size_t>()(hash_size_t_pair(p.e_.hash(), p.position_));
    }
};
}

inline std::ostream &operator<<(std::ostream &os, const DebruijnGraphCursor &p) {
  if (p.is_empty()) {
    return os << "(@)";
  } else {
    return os << "(" << p.e_ << ", " << p.position_ << ")";
  }
}

inline std::ostream &operator<<(std::ostream &os, const DebruijnComponentCursor &p) {
  if (p.is_empty()) {
    return os << "(@)";
  } else {
    return os << "(" << p.e_ << ", " << p.position_ << ")";
  }
}

std::vector<DebruijnGraphCursor> all(const debruijn_graph::ConjugateDeBruijnGraph &g);
std::vector<DebruijnComponentCursor> all(const omnigraph::GraphComponent<debruijn_graph::ConjugateDeBruijnGraph> &c);

namespace hmm {
struct Fees;
};

PathSet<DebruijnGraphCursor> find_best_path(const hmm::Fees &fees,
                                            const std::vector<DebruijnGraphCursor> &initial);
PathSet<DebruijnComponentCursor> find_best_path(const hmm::Fees &fees,
                                                const std::vector<DebruijnComponentCursor> &initial);
PathSet<ReversalGraphCursor<DebruijnGraphCursor>> find_best_path_rev(const hmm::Fees &fees,
                                                                     const std::vector<ReversalGraphCursor<DebruijnGraphCursor>> &initial);
PathSet<ReversalGraphCursor<DebruijnComponentCursor>> find_best_path_rev(const hmm::Fees &fees,
                                                                         const std::vector<ReversalGraphCursor<DebruijnComponentCursor>> &initial);
PathSet<AAGraphCursor<DebruijnComponentCursor>> find_best_path(const hmm::Fees &fees,
                                                               const std::vector<AAGraphCursor<DebruijnComponentCursor>> &initial);
