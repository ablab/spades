//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once


#include "io_support.hpp"


namespace path_extend {

using namespace debruijn_graph;


class ContigWriter {
protected:
    DECL_LOGGER("PathExtendIO")

protected:
    const Graph& g_;
    ContigConstructor<Graph> &constructor_;
    map<EdgeId, ExtendedContigIdT> ids_;
    shared_ptr<ContigNameGenerator> name_generator_;

    string ToFASTGPathFormat(const BidirectionalPath &path) const;


public:
    ContigWriter(const Graph& g,
                 ContigConstructor<Graph> &constructor,
                 const ConnectedComponentCounter &c_counter,
                 shared_ptr<ContigNameGenerator> name_generator) :
            g_(g),
            constructor_(constructor),
            ids_(),
            name_generator_(name_generator) {
        MakeContigIdMap(g_, ids_, c_counter, "NODE");
    }

    void OutputPaths(const PathContainer &paths,
                               const string &filename_base,
                               bool write_fastg = true) const;

};


class PathInfoWriter {
protected:
    DECL_LOGGER("PathExtendIO")

public:

    void WritePaths(const PathContainer &paths, const string &filename) const;
};

}
