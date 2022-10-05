//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

// This one is temporary, until we will be able to untangle IDs from the graph
#include "assembly_graph/core/graph.hpp"
#include "io/utils/id_mapper.hpp"

#include "adt/iterator_range.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace debruijn_graph {
class DeBruijnGraph;
};

namespace io {
template<typename IdType>
class IdMapper;
}

namespace gfa {
struct path;
struct segment;
struct link;



class GFAReader {
    typedef debruijn_graph::DeBruijnGraph Graph;
    typedef Graph::EdgeId EdgeId;
    typedef std::vector<std::shared_ptr<gfa::link>> Links;

  public:
    struct GFAPath {
        GFAPath(std::string n = "")
                : name(std::move(n)) {}

        std::string name;
        std::vector<EdgeId> edges;
    };
    typedef std::vector<GFAPath>::const_iterator path_const_iterator;
    typedef std::vector<GFAPath>::iterator path_iterator;

    GFAReader(const std::filesystem::path &filename)
            : filename_(filename) {}

    size_t num_edges() const { return num_edges_; }
    size_t num_links() const { return num_links_; }
    size_t num_paths() const { return paths_.size(); }

    path_const_iterator path_begin() const { return paths_.begin(); }
    path_const_iterator path_end() const { return paths_.end(); }
    adt::iterator_range<path_const_iterator> paths() const {
        return adt::make_range(path_begin(), path_end());
    }

    path_iterator path_begin() { return paths_.begin(); }
    path_iterator path_end() { return paths_.end(); }
    adt::iterator_range<path_iterator> paths() {
        return adt::make_range(path_begin(), path_end());
    }

    unsigned to_graph(debruijn_graph::DeBruijnGraph &g, io::IdMapper<std::string> *id_mapper = nullptr);

  private:
    void HandlePath(const gfa::path &record,
                    const io::IdMapper<std::string> &mapper,
                    const debruijn_graph::Graph &g);

    void HandleLink(const gfa::link &record);

    void ProcessLinks(debruijn_graph::DeBruijnGraph &g, const io::IdMapper<std::string> &mapper);

    Links link_storage_;
    std::filesystem::path filename_;
    std::vector<GFAPath> paths_;
    size_t num_edges_ = 0;
    size_t num_links_ = 0;
};

};
