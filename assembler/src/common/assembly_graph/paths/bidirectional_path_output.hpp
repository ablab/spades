//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once


#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/graph_support/contig_output.hpp"
#include "assembly_graph/components/connected_component.hpp"


namespace path_extend {

using namespace debruijn_graph;

struct IOContigStorage {
    std::string sequence_;
    BidirectionalPath* path_;

    IOContigStorage(const std::string& sequence, BidirectionalPath* path) :
    sequence_(sequence), path_(path) { }
};

struct IOContigStorageGreater
{
    bool operator()(const IOContigStorage &a, const IOContigStorage &b) const {
        if (a.sequence_.length() ==  b.sequence_.length())
            return math::gr(a.path_->Coverage(), b.path_->Coverage());
        return a.sequence_.length() > b.sequence_.length();
    }
};

//Finds common long edges in paths and joins them into
//Based on disjoint set union
class TranscriptToGeneJoiner {
private:
    const Graph &g_;
    size_t min_edge_len_; //minimal length for joining transcripts into a gene

    map<BidirectionalPath *, size_t, PathComparator> path_id_; //path ids
    std::vector<size_t> parents_; //node parents in
    std::vector<size_t> ranks_; //tree depth


    void MakeSet(size_t x);

    void JoinTrees(size_t x, size_t y);

    void Init(const PathContainer &paths);
public:
    TranscriptToGeneJoiner(const Graph &g, size_t min_edge_len): g_(g), min_edge_len_(min_edge_len) {}

    size_t FindTree(size_t x);

    size_t GetPathId(BidirectionalPath *path);

    void Construct(const PathContainer &paths);
};

class ContigWriter {
protected:
    DECL_LOGGER("PathExtendIO")

protected:
    const Graph& g_;
    ContigConstructor<Graph> &constructor_;
    size_t k_;
    map<EdgeId, ExtendedContigIdT> ids_;
    const ConnectedComponentCounter &c_counter_;
    config::pipeline_type pipeline_type_;

    string ToString(const BidirectionalPath& path) const;

    string ToFASTGString(const BidirectionalPath& path) const;


public:
    ContigWriter(const Graph& g,
                 ContigConstructor<Graph> &constructor,
                 const ConnectedComponentCounter &c_counter,
                 config::pipeline_type pipeline_type) :
            g_(g), constructor_(constructor), k_(g.k()),
            ids_(), c_counter_(c_counter),
            pipeline_type_(pipeline_type) {
        MakeContigIdMap(g_, ids_, c_counter, "NODE");
    }

    void WritePathsToFASTA(const PathContainer &paths,
                               const string &filename_base,
                               bool write_fastg = true, size_t long_edge_threshold = 300) const;

    void OutputPaths(const PathContainer& paths, const string& filename_base) const;


};


class PathInfoWriter {
protected:
    DECL_LOGGER("PathExtendIO")

public:

    void WritePaths(const PathContainer &paths, const string &filename) const;
};


class ScaffoldBreaker {
private:

    int min_gap_;

    void SplitPath(const BidirectionalPath& path, PathContainer &result) const;

public:

    ScaffoldBreaker(int min_gap): min_gap_(min_gap) {}

    void Break(const PathContainer &paths, PathContainer &result) const;
};

}
