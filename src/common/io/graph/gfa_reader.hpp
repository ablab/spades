//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
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

    typedef std::vector<std::pair<EdgeId, EdgeId>> GapLinks;
    typedef GapLinks::const_iterator jump_const_iterator;
    typedef GapLinks::iterator jump_iterator;

    size_t num_edges() const { return num_edges_; }
    size_t num_links() const { return num_links_; }
    size_t num_paths() const { return paths_.size(); }
    size_t num_gaplinks() const { return gap_links_.size(); }

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

    jump_const_iterator jump_begin() const { return gap_links_.begin(); }
    jump_const_iterator jump_end() const { return gap_links_.end(); }
    adt::iterator_range<jump_const_iterator> jumps() const {
        return adt::make_range(jump_begin(), jump_end());
    }

    jump_iterator jump_begin() { return gap_links_.begin(); }
    jump_iterator jump_end() { return gap_links_.end(); }
    adt::iterator_range<jump_iterator> jumps() {
        return adt::make_range(jump_begin(), jump_end());
    }

    unsigned to_graph(debruijn_graph::DeBruijnGraph &g, io::IdMapper<std::string> *id_mapper = nullptr);

  private:
    std::filesystem::path filename_;
    std::vector<GFAPath> paths_;
    GapLinks gap_links_;
    size_t num_edges_ = 0;
    size_t num_links_ = 0;
};

};
