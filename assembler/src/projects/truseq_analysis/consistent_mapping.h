//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace alignment_analysis {
    typedef debruijn_graph::DeBruijnGraph Graph;
    typedef Graph::VertexId VertexId;
    typedef Graph::EdgeId EdgeId;

    struct EdgeRange {
        EdgeRange(const EdgeId &first, const Range &second) : first(first), second(second) { }
        EdgeId first;
        Range second;
    };

    ostream & operator<<(ostream& os, const EdgeRange& er);

    struct ConsistentMapping {
        ConsistentMapping(const Graph &graph);

        ConsistentMapping(const Graph &graph, EdgeId e, omnigraph::MappingRange m);

        ConsistentMapping(const Graph &graph, const omnigraph::MappingPath<EdgeId> &path);

        ConsistentMapping(Graph const &graph, Range r, const vector<EdgeRange> &path);

        bool CheckConnect(EdgeId e, Range r) const;

        bool CheckConnect(const EdgeRange &er) const;

        bool CheckConnect(EdgeId e, omnigraph::MappingRange r) const;

        bool CheckConnect(const ConsistentMapping &other) const;

        bool IsEmpty() const;

        void Join(const ConsistentMapping &other);

        void Join(const ConsistentMapping &other, const vector <EdgeRange> &path);

        void ForceJoin(const ConsistentMapping &other, const vector <EdgeId> &path);

        Range const &GetInitialRange() const;

        const vector <EdgeRange> &GetMappedPath() const;

        VertexId StartVertex() const;

        VertexId EndVertex() const;

        EdgeId StartEdge() const;

        EdgeId EndEdge() const;

        const EdgeRange &Back() const;

        const EdgeRange &Front() const;

        void CutToVertex(VertexId path_start);

        Sequence CorrectSequence() const;

        size_t size() const;

        string description_;

        string CompareToReference(const string &reference_part) const;

//        void CloseEnd();
//
//        void CloseStart();

    private:
        bool CheckConnect(const EdgeRange &r1, const EdgeRange &r2) const;
        bool CheckConnect(const vector <EdgeRange> &path) const;
        vector<EdgeRange> GenerateMappingPath(const vector<EdgeId> &path) const;

        const Graph &graph_;
        Range initial_range;
        vector <EdgeRange> mapped_path;
        DECL_LOGGER("ConsistentMapping");
    };

    ostream & operator<<(ostream& os, const ConsistentMapping& cm);
}
