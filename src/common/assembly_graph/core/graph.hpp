//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "adt/iterator_range.hpp"
#include "id_storage.hpp"
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
    typedef omnigraph::impl::LinkId LinkId;
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
        LinkId conjugate_;
    };

    static constexpr unsigned LINK_ID_BIAS = 3;
    omnigraph::IdStorage<Link> lstorage_;

public:
    DeBruijnGraph(unsigned k)
            : base(k), coverage_index_(*this), lstorage_(LINK_ID_BIAS) {}

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

    using base::conjugate;

    LinkId add_link(EdgeId e1, EdgeId e2, unsigned ovl) {
        LinkId fwd = lstorage_.create(e1, e2, ovl);
        if (e1 == this->conjugate(e2)) {
            lstorage_.at(fwd.int_id()).conjugate_ = fwd;
            return fwd;
        }
        LinkId conj = lstorage_.create(this->conjugate(e2), this->conjugate(e1), ovl);
        lstorage_.at(fwd.int_id()).conjugate_ = conj;
        lstorage_.at(conj.int_id()).conjugate_ = fwd;
        return fwd;
    }

    LinkId conjugate(LinkId id) const {
        return lstorage_.at(id.int_id()).conjugate_;
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
                                     return lstorage_.at(link_id.int_id()).link.first == e;
                                   }), links.end());
    }

    void erase_links_with_outedge(VertexId v, EdgeId e) {
        auto &links = data(v).links();
        links.erase(std::remove_if(links.begin(),
                                   links.end(),
                                   [this, &e](const LinkId &link_id) {
                                       return lstorage_.at(link_id.int_id()).link.second == e;
                                   }), links.end());
    }

    void replace_outedge_in_links(VertexId v, EdgeId old_edge, EdgeId new_edge) {
        for (const LinkId &link_id : data(v).links()) {
            auto &lnk = lstorage_.at(link_id.int_id());
            if (lnk.link.second == old_edge) {
                lnk.link.second = new_edge;
                auto &conj_lnk = lstorage_.at(lnk.conjugate_.int_id());
                conj_lnk.link.first = this->conjugate(new_edge);
            }
        }
    }

    void replace_inedge_in_links(VertexId v, EdgeId old_edge, EdgeId new_edge) {
        for (const LinkId &link_id : data(v).links()) {
            auto &lnk = lstorage_.at(link_id.int_id());
            if (lnk.link.first == old_edge) {
                lnk.link.first = new_edge;
                auto &conj_lnk = lstorage_.at(lnk.conjugate_.int_id());
                conj_lnk.link.second = this->conjugate(new_edge);
            }
        }
    }

    auto links(VertexId v) const {
        return data(v).links();
    }

    const auto& link(LinkId idx) const {
        return lstorage_.at(idx.int_id());
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

    void lreserve(size_t size) { lstorage_.reserve(size); }

    size_t link_size() const { return lstorage_.size(); }

    VertexData ConjugateData(const VertexData &data) const {
        return master().conjugate(data, [this](LinkId lid) { return conjugate(lid); });
    }

    EdgeId MergePath(const std::vector<EdgeId> &path,
                     bool safe_merging = true,
                     std::vector<uint32_t> overlaps = std::vector<uint32_t>()) {
        EdgeId first_edge = path.front();
        EdgeId last_edge = path.back();
        VertexId orig_start = this->EdgeStart(first_edge);
        VertexId orig_end = this->EdgeEnd(last_edge);

        EdgeId result = base::MergePath(path, safe_merging, std::move(overlaps));

        // CorrectMergePath may conjugate-reverse the path; detect via start vertex.
        VertexId v_start = this->EdgeStart(result);
        EdgeId replaced_out = (v_start == orig_start) ? first_edge : this->conjugate(last_edge);
        if (is_complex(v_start))
            replace_outedge_in_links(v_start, replaced_out, result);

        VertexId v_end = this->EdgeEnd(result);
        EdgeId replaced_in = (v_end == orig_end) ? last_edge : this->conjugate(first_edge);
        if (is_complex(v_end))
            replace_inedge_in_links(v_end, replaced_in, result);

        return result;
    }

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
