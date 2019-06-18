//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "adt/iterator_range.hpp"
#include "observable_graph.hpp"
#include "coverage.hpp"
#include "debruijn_data.hpp"
namespace debruijn_graph {

using omnigraph::CoverageIndex;
class DeBruijnGraph: public omnigraph::ObservableGraph<DeBruijnDataMaster> {
public:
    typedef omnigraph::ObservableGraph<DeBruijnDataMaster> base;
    typedef base::DataMasterT DataMasterT;
    typedef base::VertexData VertexData;
    typedef base::EdgeData EdgeData;
    typedef base::EdgeId EdgeId;
    typedef base::VertexId VertexId;
    typedef base::VertexIt VertexIt;
    typedef DataMasterT::LinkId LinkId;
    typedef VertexIt VertexIterator;
    typedef VertexIterator iterator; // for for_each
    typedef const VertexIterator const_iterator; // for for_each
private:
    CoverageIndex<DeBruijnGraph> coverage_index_;

    struct Link {
        Link() = default;
        Link(std::pair<EdgeId, EdgeId> link, unsigned overlap)
                : link(std::move(link)), overlap(overlap) {}
        Link(EdgeId e1, EdgeId e2, unsigned overlap)
                : link{e1, e2}, overlap(overlap) {}

        Link(const Link&) = default;
        Link(Link&&) = default;

        std::pair<EdgeId, EdgeId> link;
        unsigned overlap;
    };

    typedef std::vector<Link> LinkStorage;
    LinkStorage link_storage_;

public:
    DeBruijnGraph(unsigned k)
            : base(k), coverage_index_(*this), link_storage_{} {}

    CoverageIndex<DeBruijnGraph>& coverage_index() {
        return coverage_index_;
    }

    const CoverageIndex<DeBruijnGraph>& coverage_index() const {
        return coverage_index_;
    }

    /**
     * Method returns average coverage of the edge
     */
    double coverage(EdgeId edge) const {
        return coverage_index_.coverage(edge);
    }

    uint64_t kmer_multiplicity(EdgeId edge) const {
        return coverage_index_.RawCoverage(edge);
    }

    void set_overlap(VertexId v, unsigned ovl) {
        data(v).set_overlap(ovl);
    }

    LinkId add_link(EdgeId e1, EdgeId e2, unsigned ovl) {
        link_storage_.emplace_back(e1, e2, ovl);
        return link_storage_.size() - 1;
    }

    void add_link(VertexId v, LinkId idx) {
        data(v).add_link(idx);
    }

    void add_links(VertexId v, const std::vector<LinkId> &links) {
        data(v).add_links(links);
    }

    auto move_links(VertexId v) {
        return data(v).move_links();
    }

    auto clear_links(VertexId v) {
        data(v).clear_links();
    }

    void erase_links_with_inedge(VertexId v, EdgeId e) {
        auto &links = data(v).links();
        links.erase(std::remove_if(links.begin(),
                                   links.end(),
                                   [this, &e](const LinkId &link_id) {
                                     return link_storage_[link_id].link.first == e;
                                   }), links.end());
    }

    void erase_links_with_outedge(VertexId v, EdgeId e) {
        auto &links = data(v).links();
        links.erase(std::remove_if(links.begin(),
                                   links.end(),
                                   [this, &e](const LinkId &link_id) {
                                       return link_storage_[link_id].link.second == e;
                                   }), links.end());
    }

    auto links(VertexId v) const {
        return data(v).links();
    }

    const auto& link(size_t idx) const {
        return link_storage_[idx];
    }

    bool is_complex(VertexId v) const {
        return data(v).has_complex_overlap();
    }

    size_t link_length(VertexId v, EdgeId in, EdgeId out) const {
        //todo optimize
        if (not is_complex(v)) {
            return data(v).overlap();
        }
        for (const LinkId &link_id: links(v)) {
            const Link &link = this->link(link_id);
            if (link.link.first == in and link.link.second == out) {
                return link.overlap;
            }
        }
        VERIFY_MSG(false, "Link " << in.int_id() << " -> " << out.int_id() << " was not found for vertex " << v.int_id());
    }

    void lreserve(size_t size) { link_storage_.reserve(size); }

    auto link_begin() { return link_storage_.begin(); }
    auto link_end() { return link_storage_.end(); }
    auto link_begin() const { return link_storage_.begin(); }
    auto link_end() const { return link_storage_.end(); }
    auto links() { return adt::make_range(link_begin(), link_end()); }
    auto links() const { return adt::make_range(link_begin(), link_end()); }

    size_t link_size() const {return link_storage_.size(); }

    using base::AddVertex;
    using base::AddEdge;

    VertexId AddVertex(unsigned ovl = -1U) {
        return AddVertex(VertexData(ovl == -1U ? k() : ovl));
    }

    EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
        VERIFY(nucls.size() > k());
        return AddEdge(from, to, EdgeData(nucls));
    }

    unsigned k() const {
        return master().k();
    }

    /**
     * Method returns Sequence stored in the edge
     */
    const Sequence& EdgeNucls(EdgeId edge) const {
        return this->data(edge).nucls();
    }

    const Sequence VertexNucls(VertexId v) const {
        //todo add verify on vertex nucls consistency
        if (this->OutgoingEdgeCount(v) > 0) {
            return EdgeNucls(*(this->out_begin(v))).Subseq(0, k());
        } else if (this->IncomingEdgeCount(v) > 0) {
            EdgeId inc = *(this->in_begin(v));
            size_t length = EdgeNucls(inc).size();
            return EdgeNucls(inc).Subseq(length - k(), length);
        }
        VERIFY(false);
        return Sequence();
    }

private:
    DECL_LOGGER("DeBruijnGraph")
};

typedef DeBruijnGraph ConjugateDeBruijnGraph;

typedef ConjugateDeBruijnGraph Graph;
typedef Graph::EdgeId EdgeId;
typedef Graph::VertexId VertexId;
typedef Graph::LinkId LinkId;
}
