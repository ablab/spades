//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <memory>
#include <string>

extern "C" {
    typedef struct gfa_s gfa_t;
};

namespace debruijn_graph {
class DeBruijnGraph;
};

namespace gfa {

class GFAReader {
  public:
    GFAReader();
    GFAReader(const std::string &filename);
    bool open(const std::string &filename);
    bool valid() const { return (bool)gfa_; }
    gfa_t *get() const { return gfa_.get(); }

    uint32_t num_edges() const;
    uint64_t num_links() const;

    void to_graph(debruijn_graph::DeBruijnGraph &g,
                  bool numeric_ids = true);
    
  private:
    std::unique_ptr<gfa_t, void(*)(gfa_t*)> gfa_;
};

};


