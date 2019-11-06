//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

// This one is temporary, until we will be able to untangle IDs from the graph
#include "assembly_graph/core/graph.hpp"
#include "io/utils/id_mapper.hpp"

#include "adt/iterator_range.hpp"

#include <memory>
#include <string>
#include <vector>

extern "C" {
    typedef struct gfa_s gfa_t;
};

namespace debruijn_graph {
class DeBruijnGraph;
};

namespace io {
template<typename IdType>
class IdMapper;
}

namespace gfa {

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
    typedef std::vector<GFAPath>::const_iterator path_iterator;

    GFAReader();
    GFAReader(const std::string &filename);
    bool open(const std::string &filename);
    bool valid() const { return (bool)gfa_; }
    gfa_t *get() const { return gfa_.get(); }

    uint32_t num_edges() const;
    uint64_t num_links() const;

    size_t num_paths() const { return paths_.size(); }
    path_iterator path_begin() const { return paths_.begin(); }
    path_iterator path_end() const { return paths_.end(); }
    adt::iterator_range<path_iterator> paths() const {
        return adt::make_range(path_begin(), path_end());
    }

    unsigned k() const;
    void to_graph(debruijn_graph::DeBruijnGraph &g, io::IdMapper<std::string> *id_mapper = nullptr);

  private:
    std::unique_ptr<gfa_t, void(*)(gfa_t*)> gfa_;
    std::vector<GFAPath> paths_;
};

};
