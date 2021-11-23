#include "assembly_graph/core/graph.hpp"
#include "io/graph/gfa_reader.hpp"

#include "gfa1/gfa.h"

using namespace debruijn_graph;

namespace cont_index {

class MultiplexGFAReader {
    //todo merge with GFAReader
  public:
    struct GFAPath {
      GFAPath(std::string n = "")
          : name(std::move(n)) {}

      std::string name;
      std::vector<EdgeId> edges;
    };
    typedef std::vector<GFAPath>::const_iterator path_iterator;
    typedef std::unordered_map<std::string, Sequence> OverlapStorage;
    typedef std::pair<OverlapStorage, OverlapStorage> OverlapStorages;

    MultiplexGFAReader();
    MultiplexGFAReader(const std::string &filename);
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
    OverlapStorages ConstructOverlapStorages() const;

  private:
    static void SafeInsert(OverlapStorage &storage, const std::string &name, const Sequence &seq);

    std::unique_ptr<gfa_t, void(*)(gfa_t*)> gfa_;
    std::vector<GFAPath> paths_;
};
}