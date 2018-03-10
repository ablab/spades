#include "cursor.hpp"
#include "hmmpath.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/components/graph_component.hpp"

using namespace debruijn_graph;

class DebruijnGraphCursor : public AbstractGraphCursor<DebruijnGraphCursor> {
  public:
    DebruijnGraphCursor()
            : g_{nullptr}, e_(), position_(-1) {}

    DebruijnGraphCursor(const ConjugateDeBruijnGraph *graph,
                        EdgeId e, size_t position)
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

    bool is_empty() const { return e_.get() == nullptr; }

    char letter() const { return g_->EdgeNucls(e_)[position_]; }

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

    std::vector<DebruijnGraphCursor> prev() const {
        // Case 1: edge is a tip and we're inside the terminal vertex
        if (position_ == 0) {
            // assert(pg_->ingoing_[edge_id_].size() == 0);
            return {};
        }

        // Case 2: move backwards possibly going inside the terminal vertex of a
        // tip
        if (position_ != g_->k() ||
            g_->IncomingEdgeCount(g_->EdgeStart(e_)) == 0) {
            return { DebruijnGraphCursor(g_, e_, position_ - 1) };
        }

        // Case 3: go into incoming edges
        assert(position_ == g_->k());
        assert(g_->IncomingEdgeCount(g_->EdgeStart(e_)));
        std::vector<DebruijnGraphCursor> result;
        for (EdgeId in : g_->IncomingEdges(g_->EdgeStart(e_)))
            result.emplace_back(g_, in, g_->length(in) - 1);

        return result;
    }


    std::vector<DebruijnGraphCursor> next() const {
        // Common case: we have not reached the end of the edge (in nucls)
        if (position_ + 1 < g_->length(e_) + g_->k())
            return { DebruijnGraphCursor(g_, e_, position_ + 1) };

        // Otherwise we're inside the vertex and need to go out of it
        assert(position_ + 1 == g_->length(e_) + g_->k() + 1);
        std::vector<DebruijnGraphCursor> result;
        result.reserve(4);
        for (EdgeId out : g_->OutgoingEdges(g_->EdgeEnd(e_)))
            result.emplace_back(g_, out, g_->k());  // Vertices are k-mers

        return result;
    }

  private:
    void normalize_prefix_to_suffix() {
        while (position_ < g_->k() && g_->IncomingEdgeCount(g_->EdgeStart(e_)) > 0) {
            e_ = *g_->in_begin(g_->EdgeStart(e_));
            position_ = g_->length(e_) + position_;
        }
    }

    const ConjugateDeBruijnGraph *g_;
    EdgeId e_;
    size_t position_;
};

class DebruijnComponentCursor : public AbstractGraphCursor<DebruijnComponentCursor> {
  public:
    DebruijnComponentCursor()
            : c_{nullptr}, e_(), position_(-1) {}

    DebruijnComponentCursor(const omnigraph::GraphComponent<ConjugateDeBruijnGraph> *component,
                            EdgeId e, size_t position)
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

    bool is_empty() const { return e_.get() == nullptr; }

    char letter() const { return c_->g().EdgeNucls(e_)[position_]; }

    bool is_convergent() const {
        const ConjugateDeBruijnGraph &g = c_->g();
        return position_ == g.k() &&
                (g.IncomingEdgeCount(g.EdgeStart(e_)) > 1 && !c_->IsBorder(g.EdgeStart(e_)));
        // return prev().size() > 1;
    }

    bool is_divergent() const {
        const ConjugateDeBruijnGraph &g = c_->g();
        return position_ + 1 == g.length(e_) + g.k() &&
                (g.OutgoingEdgeCount(g.EdgeEnd(e_)) > 1 && !c_->IsBorder(g.EdgeEnd(e_)));
        // return next().size() > 1;
    }

    std::vector<DebruijnComponentCursor> prev() const {
        const ConjugateDeBruijnGraph &g = c_->g();

        // Case 1: edge is a tip and we're inside the terminal vertex
        if (position_ == 0) {
            // assert(pc_->g().ingoing_[edge_id_].size() == 0);
            return {};
        }

        // Case 2: move backwards possibly going inside the terminal vertex of a
        // tip
        if (position_ != g.k() ||
            (g.IncomingEdgeCount(g.EdgeStart(e_)) == 0 ||
             c_->IsBorder(g.EdgeStart(e_)))) {
            return { DebruijnComponentCursor(c_, e_, position_ - 1) };
        }

        // Case 3: go into incoming edges
        assert(position_ == g.k());
        assert(g.IncomingEdgeCount(g.EdgeStart(e_)) && !c_->IsBorder(g.EdgeStart(e_)));
        std::vector<DebruijnComponentCursor> result;
        for (EdgeId in : g.IncomingEdges(g.EdgeStart(e_)))
            result.emplace_back(c_, in, g.length(in) - 1);

        return result;
    }


    std::vector<DebruijnComponentCursor> next() const {
        const ConjugateDeBruijnGraph &g = c_->g();

        // Common case: we have not reached the end of the edge (in nucls)
        if (position_ + 1 < g.length(e_) + g.k())
            return { DebruijnComponentCursor(c_, e_, position_ + 1) };

        // If we're in the border vertex, then no extension is possible
        if (c_->IsBorder(g.EdgeEnd(e_)))
            return { };

        // Otherwise we're inside the vertex and need to go out of it
        assert(position_ + 1 == g.length(e_) + g.k() + 1);
        std::vector<DebruijnComponentCursor> result;
        result.reserve(4);
        for (EdgeId out : g.OutgoingEdges(g.EdgeEnd(e_)))
            result.emplace_back(c_, out, g.k());  // Vertices are k-mers

        return result;
    }

  private:
    void normalize_prefix_to_suffix() {
        const ConjugateDeBruijnGraph &g = c_->g();
        while (position_ < g.k() &&
               (g.IncomingEdgeCount(g.EdgeStart(e_)) > 0 && !c_->IsBorder(g.EdgeStart(e_)))) {
            e_ = *g.in_begin(g.EdgeStart(e_));
            position_ = g.length(e_) + position_;
        }
    }

    const omnigraph::GraphComponent<ConjugateDeBruijnGraph> *c_;
    EdgeId e_;
    size_t position_;
};
