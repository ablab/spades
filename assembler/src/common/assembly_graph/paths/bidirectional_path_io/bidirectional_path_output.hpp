//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/utils/edge_namer.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/graph/gfa_writer.hpp"
#include "io/graph/fastg_writer.hpp"
#include "io_support.hpp"

#include <unordered_set>

namespace path_extend {

class PathWriter {
    typedef debruijn_graph::DeBruijnGraph Graph;

  public:
    PathWriter(const Graph &graph)
            : graph_(graph),
              short_namer_(graph) {}

    std::string ToPathString(const BidirectionalPath &path) const;

  private:
    const Graph &graph_;
    io::CanonicalEdgeHelper<Graph> short_namer_;
};

class FastgPathWriter : public io::FastgWriter {
public:
    using io::FastgWriter::FastgWriter;

    FastgPathWriter(const Graph &graph,
                    const std::filesystem::path &fn,
                    io::EdgeNamingF<Graph> edge_naming_f = io::BasicNamingF<Graph>())
            :  io::FastgWriter(graph, fn, edge_naming_f),
               path_writer_(graph)
    {}

    void WritePaths(const ScaffoldStorage &scaffold_storage, const std::filesystem::path &fn) const;

  private:
    PathWriter path_writer_;
};


template<class Base>
class GFAPathWriterBase : public Base {
    void WritePath(const std::string &name,
                   const std::vector<std::string> &edge_strs,
                   const std::string &flags) {
        this->os_ << "P" << "\t" ;
        this->os_ << name << "\t";
        bool first = true;
        for (const auto& e : edge_strs) {
            this->os_ << (first ? "" : "," ) << e;
            first = false;
        }
        this->os_ << "\t*";
        if (flags.length())
            this->os_ << '\t' << flags;
        this->os_ << '\n';
    }

    void WritePath(const std::string &name, size_t segment_id,
                   const std::vector<std::string> &edge_strs,
                   const std::string &flags = "");

    void WritePath(const std::string &name,
                   const std::vector<std::string> &edge_strs,
                   const std::string &flags = "");

    void WritePaths11(const std::vector<EdgeId> &edges,
                      const std::string &name,
                      const std::string &flags = "");

    void WritePaths11(const ScaffoldStorage &scaffold_storage);

    void WritePaths12(const std::vector<EdgeId> &edges,
                      const std::string &name,
                      const std::string &flags = "");

    void WritePaths12(const ScaffoldStorage &scaffold_storage);

    using JumpLinks = std::unordered_set<std::pair<EdgeId, EdgeId>>;
    void WriteJumpLinks(const JumpLinks &links);

public:
    enum class Version {
        GFAv11, // Split gapped paths into segments
        GFAv12  // Add jump links
    };

    GFAPathWriter(const Graph &graph, std::ostream &os,
                  io::EdgeNamingF<Graph> naming_f = io::IdNamingF<Graph>(),
                  Version version = Version::GFAv12);

    void WritePaths(const std::vector<EdgeId> &edges,
                    const std::string &name,
                    const std::string &flags = "");

    void WritePaths(const ScaffoldStorage &scaffold_storage);

    void WritePaths(const gfa::GFAReader::GFAPath &path,
                    const std::string &flags = "");

private:
    Version version_;
};

using GFAPathWriter = GFAPathWriterBase<gfa::GFAWriter>;
using GFAPathComponentWriter = GFAPathWriterBase<gfa::GFAComponentWriter>;


typedef std::function<void (const ScaffoldStorage&)> PathsWriterT;

class ContigWriter {
    const Graph& g_;
    std::shared_ptr<ContigNameGenerator> name_generator_;

public:
    static void WriteScaffolds(const ScaffoldStorage &scaffold_storage, const std::filesystem::path &fn) {
        io::OFastaReadStream oss(fn);
        std::ofstream os_fastg;

        for (const auto& scaffold_info : scaffold_storage) {
            TRACE("Scaffold " << scaffold_info.name << " originates from path " << scaffold_info.path->str());
            oss << io::SingleRead(scaffold_info.name, scaffold_info.sequence);
        }
    }

    static PathsWriterT BasicFastaWriter(const std::filesystem::path &fn) {
        return [=](const ScaffoldStorage& scaffold_storage) {
            WriteScaffolds(scaffold_storage, fn);
        };
    }

    ContigWriter(const Graph& g,
                 std::shared_ptr<ContigNameGenerator> name_generator) :
            g_(g),
            name_generator_(name_generator) {
    }

    void OutputPaths(const PathContainer &paths, const std::vector<PathsWriterT>& writers) const;

    void OutputPaths(const PathContainer &paths, PathsWriterT writer) const {
        OutputPaths(paths, std::vector<PathsWriterT>{writer});
    }

    void OutputPaths(const PathContainer &paths, const std::filesystem::path &fn) const {
        OutputPaths(paths, BasicFastaWriter(fn));
    }

private:
    DECL_LOGGER("ContigWriter")
};

}
