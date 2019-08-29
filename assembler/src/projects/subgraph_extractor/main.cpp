//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "subgraph_extraction.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "pipeline/config_struct.hpp"

#include "toolchain/edge_label_helper.hpp"
#include "toolchain/utils.hpp"

#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include <clipp/clipp.h>
#include <string>
#include <set>
#include <numeric>
#include <sys/types.h>

using namespace debruijn_graph;

struct gcfg {
    gcfg()
        : k(0), //save_gfa(false),
          nthreads(omp_get_max_threads() / 2 + 1)
    {}

    unsigned k;
    std::string graph;
    std::string tmpdir;
    std::string outdir;
    std::string genes_desc;
    std::string genes_seq;
    std::string cds_len_fn;
//    bool save_gfa;
    unsigned nthreads;
//    output_type mode;
};

static void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
//      cfg.graph << value("graph. In GFA (ending with .gfa) or prefix to SPAdes graph pack"),
      //cfg.outfile << value("output filename/prefix (in case of --spades-gp)"),
//      option("--gfa").set(cfg.save_gfa) % "produce GFA output",
      (required("-o") & value("dir", cfg.outdir)) % "outpur directory to use for GFA files",
      one_of((option("-part-desc") & value("file", cfg.genes_desc)) % "file with partial genes description (.gff)",
             (option("-part-seq") & value("file", cfg.genes_seq)) % "file with partial genes sequences (.fasta)"),
      (required("-graph") & value("graph", cfg.graph)) % "In GFA (ending with .gfa) or prefix to SPAdes graph pack",
      (required("-cds-len-est") & value("file", cfg.cds_len_fn)) % "file with cds length estimamtes",
      (required("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use (default: max_threads / 2)",
      (option("-tmpdir") & value("dir", cfg.tmpdir)) % "scratch directory to use (default: <outdir>/tmp)"
//      one_of(option("--unitigs").set(cfg.mode, output_type::unitigs) % "produce unitigs (default)",
//             option("--fastg").set(cfg.mode, output_type::fastg) % "produce graph in FASTG format",
//             option("--gfa").set(cfg.mode, output_type::gfa) % "produce graph in GFA1 format",
//             option("--spades").set(cfg.mode, output_type::spades) % "produce graph in SPAdes internal format",
//             option("--spades-gp").set(cfg.mode, output_type::spades_pack) % "produce graph pack in SPAdes internal format "
//                                                                        "(recommended if bulges are removed to improve further read mapping)")
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }
}


//TODO simplify
using CodonSet = std::vector<Sequence>;
CodonSet STOP_CODONS = {Sequence("TAG"), Sequence("TAA"), Sequence("TGA")};

static CodonSet RC(const CodonSet &codons) {
    CodonSet rc(codons.size());
    std::transform(codons.begin(), codons.end(), rc.begin(),
                   [](const Sequence &s) {return !s;});
    return rc;
}

CodonSet RC_STOP_CODONS = RC(STOP_CODONS);;

//EdgeId + offset pair
using GraphPos = std::pair<EdgeId, size_t>;

template<typename T>
void hash_combine(std::size_t &seed, T const &key) {
    std::hash<T> hasher;
    seed ^= hasher(key) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
template<>
struct hash<GraphPos> {
    size_t operator()(const GraphPos &gpos) const {
        std::size_t seed(0);
        hash_combine(seed, gpos.first);
        hash_combine(seed, gpos.second);
        return seed;
    }
};
}

struct FramedPos {
    EdgeId e;
    size_t offset;
    //"state" of the frame: 0 means that the codon is complete
    unsigned frame;

    FramedPos(EdgeId e_, size_t offset_, unsigned frame_): e(e_), offset(offset_), frame(frame_) {}

    FramedPos(EdgeId e_, size_t offset_): FramedPos(e_, offset_, 0) {}

    FramedPos(): e(EdgeId()), offset(size_t(-1)) {}

    bool good_frame() const {
        return frame == 0;
    }

    bool last(const Graph &g) const {
        VERIFY(offset < g.length(e));
        return offset == g.length(e) - 1;
    }

    FramedPos next() const {
        return FramedPos(e, offset + 1, (frame + 1) % 3);
    }

    FramedPos next(EdgeId neighbour) const {
        return FramedPos(neighbour, 0, (frame + 1) % 3);
    }

    bool operator!=(const FramedPos &other) const {
        return e != other.e || offset != other.offset || frame != other.frame;
    }

    bool operator==(const FramedPos &other) const {
        return !(*this != other);
    }

    size_t hash() const {
        std::size_t seed(0);
        hash_combine(seed, e);
        hash_combine(seed, offset);
        hash_combine(seed, frame);
        return seed;
    }
};

namespace std {
template<>
struct hash<FramedPos> {
    std::size_t operator()(const FramedPos &fp) const {
        return fp.hash();
    }
};
}

using EdgePath = omnigraph::Path<EdgeId>;

class CodonFinder {
    const Graph &g_;
    GraphPos init_pos_;
    const CodonSet terminators_;
    std::queue<std::pair<FramedPos, FramedPos>> queue_;
    std::unordered_map<FramedPos, FramedPos> prev_;

    //Last three nucleotides of a k+1-mer on position pos in the graph
    //TODO optimize ?
    Sequence Codon(GraphPos pos) const {
        VERIFY(pos.second < g_.length(pos.first));
        return g_.EdgeNucls(pos.first).Subseq(pos.second + g_.k() - 2, pos.second + g_.k() + 1);
    }

    bool CodonInSet(GraphPos pos, const CodonSet& codons) const {
        return std::find(codons.begin(), codons.end(), Codon(pos)) != codons.end();
    }

    void NextToQueue(FramedPos fpos) {
        if (fpos.last(g_)) {
            for (EdgeId e : g_.OutgoingEdges(g_.EdgeEnd(fpos.e))) {
                queue_.push(std::make_pair(fpos.next(e), fpos));
            }
        } else {
            queue_.push(std::make_pair(fpos.next(), fpos));
        }
    }

    bool Terminate(FramedPos fpos) const {
        return fpos.good_frame() && CodonInSet(GraphPos(fpos.e, fpos.offset), terminators_);
    }

public:
    CodonFinder(const Graph &g, GraphPos init_pos, const CodonSet& terminators):
            g_(g), init_pos_(init_pos), terminators_(terminators) {}

    std::vector<EdgePath> Go() {
        FramedPos init(init_pos_.first, init_pos_.second);
        queue_.push(std::make_pair(init, FramedPos()));

        std::vector<FramedPos> terminated;

        while (!queue_.empty()) {
            FramedPos fpos;
            FramedPos prev;
            std::tie(fpos, prev) = queue_.front();
            queue_.pop();
            if (!prev_.count(fpos)) {
                prev_[fpos] = prev;
                if (Terminate(fpos)) {
                    DEBUG("Terminate graph pos: " << g_.str(fpos.e) << " " << fpos.offset);
                    DEBUG("Codon start coord " << fpos.offset + g_.k() - 1);
                    DEBUG("Codon " << Codon(GraphPos(fpos.e, fpos.offset)));
                    terminated.push_back(fpos);
                } else {
                    NextToQueue(fpos);
                }
            }
        }

        std::vector<EdgePath> paths;

        for (FramedPos t : terminated) {
            std::vector<EdgeId> reverse_path;
            reverse_path.push_back(t.e);
            if (t == init) {
                paths.push_back(EdgePath(std::vector<EdgeId>(reverse_path.rbegin(), reverse_path.rend()),
                                     init.offset + 1, init.offset + 1));
            } else {
                FramedPos fpos = t;
                while (true) {
                    VERIFY(prev_.count(fpos));
                    auto prev = prev_[fpos];
                    if (prev == init) {
                        break;
                    }
                    if (prev.e != fpos.e || fpos.offset != prev.offset + 1) {
                        reverse_path.push_back(prev.e);
                    }
                    fpos = prev;
                }
                paths.push_back(EdgePath(std::vector<EdgeId>(reverse_path.rbegin(), reverse_path.rend()),
                                     fpos.offset, t.offset + 1));
            }
        }

        return paths;
    }

    //map from graph position to its corresponding path length
    //Note that those are "exclusive" coordinates of the ends of the path
    std::unordered_map<GraphPos, size_t> Terminates(const std::vector<EdgePath> &paths) const {
        std::unordered_map<GraphPos, size_t> pos_len;
        for (const auto &p : paths) {
            pos_len[GraphPos(p.sequence().back(), p.end_pos())] = PathLength(g_, p);
        }
        return pos_len;
    };

    std::set<EdgeId> Edges(const std::vector<EdgePath> &paths) const {
        std::set<EdgeId> answer;
        for (const auto &p : paths) {
            utils::insert_all(answer, p.sequence());
        }
        return answer;
    }

};

static size_t RoundedProduct(size_t l, double coeff) {
    return size_t(math::round(coeff * double(l)));
}

template<class Graph>
class EdgeTrackingCallback: public omnigraph::PathProcessor<Graph>::Callback {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

public:
    EdgeTrackingCallback(const Graph& g) :
            g_(g) {
    }

    void HandleReversedPath(const std::vector<EdgeId>& path) override {
        utils::insert_all(edges_, path);
    }

    const std::set<EdgeId> &edges() const {
        return edges_;
    }

private:
    const Graph& g_;
    std::set<EdgeId> edges_;
};

static std::string PrintEdgePath(const Graph &g, const EdgePath &path) {
    std::stringstream ss;
    ss << "[" << omnigraph::PrintPath(g, path.sequence()) << "], ";
    ss << "start: " << path.start_pos() << ", end: " << path.end_pos();
    return ss.str();
}

class MinDistRelevantComponentFinder {
    //TODO use throughout file
    using DistInfo = std::unordered_map<GraphPos, size_t>;
    using BaseDistF = std::function<size_t (GraphPos, GraphPos)>;
    const Graph &g_;
    const double max_len_coeff_;

    size_t MaxDist(const DistInfo &v_ds) const {
        size_t max = 0;
        for (const auto &v_d : v_ds) {
            if (v_d.second > max)
                max = v_d.second;
        }
        return max;
    }

    //returns min distance to exit among the appropriate paths or max value if none
    size_t Check(const DistInfo &s_dist, const DistInfo &e_dist, BaseDistF base_dist_f) const {
        size_t min_e_dist = std::numeric_limits<size_t>::max();
        //todo seems like can be simplified
        for (const auto s_d : s_dist) {
            for (const auto e_d : e_dist) {
                if (base_dist_f(s_d.first, e_d.first) >= s_d.second + e_d.second) {
                    if (e_d.second < min_e_dist) {
                        min_e_dist = e_d.second;
                    }
                }
            }
        }
        return min_e_dist;
    }

    bool CheckConnectedTo(VertexId v, const std::set<VertexId> &vertices) const {
        for (EdgeId e : g_.IncomingEdges(v))
            if (vertices.count(g_.EdgeStart(e)))
                return true;
        return false;
    }

public:
    //TODO max_len_frac used twice with slightly different meaning
    MinDistRelevantComponentFinder(const Graph &g,
                                   double max_len_frac = 1.5) :
            g_(g), max_len_coeff_(max_len_frac) {}

    //TODO check how base_len is defined
    GraphComponent<Graph> RelevantComponent(size_t cds_len_est,
                                            size_t base_len,
                                            const DistInfo &starts,
                                            const DistInfo &ends) const {
        INFO("Extracting relevant component");
        if (starts.size() * ends.size() == 0) {
            ERROR("Set of stop codons to the left/right was empty");
            VERIFY(false);
        }
        std::unordered_map<VertexId, DistInfo> s_dists;
        std::unordered_map<VertexId, DistInfo> e_dists;
        const size_t max_cds_len = RoundedProduct(cds_len_est, max_len_coeff_);

        //TODO check Dijkstra return codes
        size_t global_bound = RoundedProduct(std::max(MaxDist(starts), MaxDist(ends)) + base_len, max_len_coeff_);
        DEBUG("Max starts dist " << MaxDist(starts));
        DEBUG("Max ends dist " << MaxDist(ends));
        DEBUG("Base len " << base_len);
        DEBUG("Global length bound " << global_bound);

        //TODO can optimize! First process direction with fewer vertices, then block unreached vertices
        //TODO we can ignore light edges
        const size_t MAX_VERTEX_BOUND = 10000;
        auto fwd_dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, global_bound,
                                                                         MAX_VERTEX_BOUND);
        DEBUG("Launching forward Dijkstras from ends of " << starts.size() << " edges");
        //TODO can be slightly optimized for several starts on same edge
        for (const auto &p : utils::key_set(starts)) {
            fwd_dijkstra.Run(g_.EdgeEnd(p.first));
            //if (fwd_dijkstra.VertexLimitExceeded()) {
            //    WARN("Vertex limit was exceeded while launching Dijkstra "
            //            "from end of edge " << g_.str(p.first));
            //    return GraphComponent<Graph>(g_);
            //}
            for (VertexId v : fwd_dijkstra.ReachedVertices()) {
                s_dists[v][p] = fwd_dijkstra.GetDistance(v) + g_.length(p.first) - p.second;
            }
        }

        if (fwd_dijkstra.VertexLimitExceeded()) {
            DEBUG("Vertex limit was exceeded while launching forward Dijkstras");
        }

        DEBUG("Launching backward Dijkstras from starts of " << ends.size() << " edges");
        //TODO deduplicate?!!
        auto bwd_dijkstra = DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(g_, global_bound,
                                                                                 MAX_VERTEX_BOUND);
        for (const auto &p : utils::key_set(ends)) {
            bwd_dijkstra.Run(g_.EdgeStart(p.first));
            //if (bwd_dijkstra.VertexLimitExceeded()) {
            //    WARN("Vertex limit was exceeded while launching Dijkstra "
            //            "from start of edge " << g_.str(p.first));
            //    return GraphComponent<Graph>(g_);
            //}
            for (VertexId v : bwd_dijkstra.ReachedVertices()) {
                e_dists[v][p] = bwd_dijkstra.GetDistance(v) + p.second;
            }
        }

        if (bwd_dijkstra.VertexLimitExceeded()) {
            DEBUG("Vertex limit was exceeded while launching backward Dijkstras");
        }

        DEBUG("Intersecting results");
        //TODO correct usage of max_len_coeff to nucleote length
        BaseDistF base_dist_f =
                [this, &starts, &ends, base_len] (GraphPos s, GraphPos e) {
                    return RoundedProduct(utils::get(starts, s) + utils::get(ends, e) + base_len, max_len_coeff_);
        };

        //vertices with minimal distance to any end among valid paths
        std::vector<std::pair<size_t, VertexId>> dist_vertices;

        for (const auto &v_di: s_dists) {
            VertexId v = v_di.first;
            const auto &s_dist = v_di.second;

            auto it = e_dists.find(v);
            if (it == e_dists.end())
                continue;

            const auto &e_dist = it->second;
            dist_vertices.push_back(std::make_pair(Check(s_dist, e_dist, base_dist_f), v));
        }

        DEBUG("Creating component");
        std::sort(dist_vertices.begin(), dist_vertices.end());
        std::set<VertexId> within_cds_limit;
        std::set<VertexId> component_vertices;

        //FIXME what if corresponding stop codon is too far?
        //Probably need to check stop codon distance first!

        //TODO add finding potential start codons?
        for (auto d_v : dist_vertices) {
            size_t min_end_dist = d_v.first;
            VertexId v = d_v.second;
            if (min_end_dist < max_cds_len) {
                within_cds_limit.insert(v);
                component_vertices.insert(v);
            } else if (min_end_dist < std::numeric_limits<size_t>::max()) {
                if (std::any_of(g_.out_begin(v), g_.out_end(v), [&](EdgeId e) {
                    return within_cds_limit.count(g_.EdgeEnd(e));
                })) {
                    component_vertices.insert(v);
                }
            }
        }

        //Adding edges containing stop codons forward (if length check passed)
        for (const auto &p_d : ends) {
            EdgeId e = p_d.first.first;
            if (base_len + p_d.second < max_cds_len) {
                if (!CheckConnectedTo(g_.EdgeStart(e), within_cds_limit))
                    WARN("Edge " << g_.str(e) << " being added as containing stop codon "
                            "is disconnected from inner component vertices");

                component_vertices.insert(g_.EdgeStart(e));
                component_vertices.insert(g_.EdgeEnd(e));
            } else {
                WARN("Edge " << g_.str(e) << " containing stop codon "
                        "on position " << p_d.first.second << " was not added to component");
                INFO("Max cds len " << max_cds_len << " ; base_len " << base_len << " ; dist to stop " << p_d.second);
            }
        }

        DEBUG("Relevant component extracted");
        return GraphComponent<Graph>::FromVertices(g_, component_vertices.begin(), component_vertices.end());
    }

};

class PathFindingRelevantComponentFinder {
    const Graph &g_;
    const double min_len_frac_;
    const double max_len_frac_;

public:
    PathFindingRelevantComponentFinder(const Graph &g, double min_len_frac = 0.5,
                         double max_len_frac = 1.5) :
            g_(g), min_len_frac_(min_len_frac), max_len_frac_(max_len_frac) {}

    GraphComponent<Graph> RelevantComponent(size_t base_len,
                                            const std::unordered_map<VertexId, size_t> &starts,
                                            const std::unordered_map<VertexId, size_t> &ends) const {
        EdgeTrackingCallback<Graph> callback(g_);
        for (auto s_v_d : starts) {
            VertexId start = s_v_d.first;
            DEBUG("Searching paths for source " << g_.str(start));
            size_t max_length = 0;
            for (auto e_v_d : ends) {
                size_t total_length = s_v_d.second + e_v_d.second + base_len;
                if (total_length > max_length)
                    max_length = total_length;
            }
            max_length = size_t(math::round(max_len_frac_ * double(max_length)));
            DEBUG("Creating path processor for max path length " << max_length);
            //todo increase dijkstra limits
            PathProcessor<Graph> processor(g_, start, max_length);
            for (auto e_v_d : ends) {
                VertexId end = e_v_d.first;
                DEBUG("Finding paths to sink " << g_.str(end));
                size_t total_length = s_v_d.second + e_v_d.second + base_len;
                processor.Process(end,
                                  size_t(math::round(min_len_frac_ * double(total_length))),
                                  size_t(math::round(max_len_frac_ * double(total_length))),
                                  callback
                );
            }
        }

        const auto &edges = callback.edges();
        return GraphComponent<Graph>::FromEdges(g_, edges.begin(), edges.end());
    }

};

class PartialGeneProcessor {
    const Graph &g_;
//    const PathFindingRelevantComponentFinder rel_comp_finder_;
    const MinDistRelevantComponentFinder rel_comp_finder_;
    const ReadPathFinder<Graph> path_finder_;
    const io::EdgeNamingF<Graph> edge_naming_f_;

    Sequence PathSequence(const EdgePath &p) const {
        return debruijn_graph::PathSequence(g_, p);
    }

    EdgePath RCEdgePath(const EdgePath &p) const {
        const auto &edges = p.sequence();
        std::vector<EdgeId> rc_edges(edges.size());
        std::transform(edges.rbegin(), edges.rend(), rc_edges.begin(), [this](EdgeId e) {return g_.conjugate(e);});
        return EdgePath(rc_edges, g_.length(edges.back()) - p.end_pos(), g_.length(edges.front()) - p.start_pos());
    }

    std::vector<EdgePath> RCEdgePaths(const std::vector<EdgePath> &ps) const {
        std::vector<EdgePath> paths(ps.size());
        std::transform(ps.begin(), ps.end(), paths.begin(), [this](const EdgePath &p) {return RCEdgePath(p);});
        return paths;
    }

    GraphPos RCPos(GraphPos p) const {
        return std::make_pair(g_.conjugate(p.first), g_.length(p.first) - p.second);
    }

    std::unordered_map<GraphPos, size_t> RCTerminates(const std::unordered_map<GraphPos, size_t> &terminate_dists) const {
        std::unordered_map<GraphPos, size_t> ans;
        for (auto v_d : terminate_dists) {
            ans[RCPos(v_d.first)] = v_d.second;
        }
        return ans;
    };

    //potential stop codons of the gene will be added to stop_codon_poss
    GraphComponent<Graph> ProcessPartialGenePath(const EdgePath &gene_path,
                                                 size_t cds_len_est,
                                                 std::set<GraphPos> &stop_codon_poss) {
        if (gene_path.size() == 0) {
            WARN("Partial gene path was empty");
            return GraphComponent<Graph>(g_);
        }
        DEBUG("Processing partial gene path " << PrintEdgePath(g_, gene_path));
        //DEBUG("Nucls: " << PathSequence(gene_path));

        INFO("Searching for stop codon forward");
        EdgeId last_e = gene_path.sequence().back();
        CodonFinder fwd_finder(g_, GraphPos(last_e, gene_path.end_pos() - 1), STOP_CODONS);
        auto paths_fwd = fwd_finder.Go();
        INFO("Codons forward found: " << paths_fwd.size());
        auto ends = fwd_finder.Terminates(paths_fwd);

        INFO("Added to set of potential stop codons for the gene");
        for (const auto &gp_d : ends) {
            //NB: Inserted number gives the (k+1-mer) coordinate beyond the
            //    k+1-mer whose sequence ends with stop codon
            stop_codon_poss.insert(gp_d.first);
        }

        TRACE("Paths forward:");
        for (const auto &p : paths_fwd) {
            TRACE("Path: " << PrintEdgePath(g_, p));
            Sequence s = PathSequence(p).Subseq(g_.k());
            TRACE("Sequence: len=" << s.size() << "; " << s);
        }

        INFO("Searching for stop codon backward");
        EdgeId first_e = gene_path.sequence().front();
        CodonFinder bwd_finder(g_, GraphPos(g_.conjugate(first_e), g_.length(first_e) - gene_path.start_pos() - 1),
                               RC_STOP_CODONS);
        auto paths_bwd = bwd_finder.Go();
        INFO("Codons backward found: " << paths_bwd.size());
        auto starts = RCTerminates(bwd_finder.Terminates(paths_bwd));

        TRACE("Paths backward (reversed): ");
        for (const auto &p : RCEdgePaths(paths_bwd)) {
            TRACE("Path: " << PrintEdgePath(g_, p));
            Sequence s = PathSequence(p);
            TRACE("Sequence: len=" << s.size() - g_.k() << "; " << s.Subseq(0, s.size() - g_.k()));
        }

        INFO("Graph sinks to consider: " << ends.size());
        //TODO optimize debug logging
        for (const auto &gpos_d: ends) {
            EdgeId e = gpos_d.first.first;
            DEBUG("Edge " << edge_naming_f_(g_, e) << " pos " << gpos_d.first.second);
        }
        //TODO optimize debug logging
        INFO("Graph sources to consider: " << starts.size());
        for (const auto &gpos_d: starts) {
            EdgeId e = gpos_d.first.first;
            DEBUG("Edge " << edge_naming_f_(g_, e) << " pos " << gpos_d.first.second);
        }

        return (starts.size() * ends.size() == 0) ? GraphComponent<Graph>(g_) :
                    rel_comp_finder_.RelevantComponent(cds_len_est, PathLength(g_, gene_path), starts, ends);
    }

//    EdgePath CropRight(const EdgePath &path, size_t min_len) const {
//        if (path.size() == 0)
//            return EdgePath();
//
//        const size_t start_pos = path.start_pos();
//        std::vector<EdgeId> edges = {path.sequence().front()};
//        size_t total_len = g_.length(path.sequence().front()) - start_pos;
//        size_t i = 1;
//
//        while (i < path.size() && total_len < min_len) {
//            EdgeId e = path[i];
//            edges.push_back(e);
//            total_len += g_.length(e);
//            ++i;
//        }
//
//        size_t end_pos;
//        if (i == path.size())
//            end_pos = path.end_pos();
//        else
//            end_pos = g_.length(path[i - 1]);
//
//        return EdgePath(edges, start_pos, end_pos);
//    }
//
//    EdgePath CropLeft(const EdgePath &path, size_t min_len) const {
//        return RCEdgePath(CropRight(RCEdgePath(path), min_len));
//    }

    template<class Mapper>
    EdgePath ExtractPath(const std::string &s, const Mapper &mapper) const {
        io::SingleRead r("empty_name", s);
        //TODO output more info about particular prediction mapping!
        EdgePath path = path_finder_.FindDetailedReadPath(mapper->MapRead(r));
        if (!CheckContiguous(g_, path.sequence())) {
//            TRACE("Path for one of predictions for " << gene_id << " not contiguous");
            TRACE("Path not contiguous");
        }
        if (PathLength(g_, path) + g_.k() < RoundedProduct(r.size(), 0.7)) {
//            WARN("Path for one of predictions for " << gene_id << " was much shorter than query");
            WARN("Path much shorter than query");
            return EdgePath();
        } else {
            return path;
        }
    }

public:
    PartialGeneProcessor(const Graph &g,
                         const io::EdgeNamingF<Graph> &edge_naming_f,
                         //double min_len_frac = 0.5,
                         double max_len_coeff = 1.5) :
            g_(g),
            rel_comp_finder_(g_, /*min_len_frac, */max_len_coeff),
            path_finder_(g, false, /*max gap length*/100),
            edge_naming_f_(edge_naming_f) {}

    //TODO configure + might need to use query length instead of path length
    //TODO think about interface
    template<class Mapper>
    GraphComponent<Graph> ProcessPartialGene(const std::string &s,
                                             const Mapper &mapper,
                                             size_t cds_length_estimate,
                                             std::set<GraphPos> &stop_codon_poss,
                                             size_t min_len_to_explore = 200,
                                             double frac_to_explore = 0.5) {
        VERIFY_MSG(s.size() % 3 == 0, "Size of partial CDS prediction not divisible by 3");
        VERIFY(math::le(frac_to_explore, 1.));
        INFO("Original query length " << s.size());
        INFO("CDS length estimate " << cds_length_estimate);
        if (s.size() < min_len_to_explore) {
            WARN("Size of partial CDS prediction shorter than " << min_len_to_explore);
            min_len_to_explore = s.size();
        }
        size_t crop_len = std::max(min_len_to_explore, RoundedProduct(s.size(), frac_to_explore));
        //aligning on frame
        crop_len = crop_len / 3 * 3;
        INFO("Will explore left/right parts of length " << crop_len);

        std::set<EdgeId> edges;

        INFO("Processing rightmost part of the query of length " << crop_len);
        //s.Last(crop_len) & s.First(crop_len)
        VERIFY(crop_len <= s.size());
        utils::insert_all(edges, ProcessPartialGenePath(ExtractPath(s.substr(s.size() - crop_len), mapper),
                    cds_length_estimate, stop_codon_poss).edges());
        INFO("Processing leftmost part of the query of length " << crop_len);
        utils::insert_all(edges, ProcessPartialGenePath(ExtractPath(s.substr(0, crop_len), mapper),
                    cds_length_estimate, stop_codon_poss).edges());

        if (edges.size() == 0) {
            INFO("Couldn't gather any reasonable component");
            return GraphComponent<Graph>(g_);
        }

        INFO("'Closing' gathered component");
        subgraph_extraction::ComponentExpander expander(g_);
        return expander.Expand(GraphComponent<Graph>::FromEdges(g_, edges.begin(), edges.end()));
    }
};

struct PartialGeneInfo {
    size_t unitig_id;
    Range r;
    std::string gene_id;
    bool strand;
};

inline std::istream& operator>>(std::istream& is, PartialGeneInfo &pgi) {
    is >> pgi.unitig_id;
    is >> pgi.r.start_pos;
    is >> pgi.r.end_pos;
    is >> pgi.gene_id;
    std::string strand_symbol;
    is >> strand_symbol;
    if (strand_symbol == "+")
        pgi.strand = true;
    else if (strand_symbol == "-")
        pgi.strand = false;
    else
        VERIFY_MSG(false, "Unsupported strand symbol");
    return is;
}

inline std::ostream& operator<<(std::ostream& os, const PartialGeneInfo &pgi) {
    os << "Unitig id: " << pgi.unitig_id << "; ";
    os << "Range: [" << pgi.r.start_pos << ", " << pgi.r.end_pos << "]; ";
    os << "Gene id:" << pgi.gene_id << "; ";
    os << "Strand: " << (pgi.strand ? "+" : "-");
    return os;
}

static Range ConjugateRange(const Graph &g, EdgeId e, const Range &r) {
    return Range(g.length(e) - r.end_pos, g.length(e) - r.start_pos);
}

using GeneInitSeq = std::multimap<std::string, std::string>;

static GeneInitSeq
PredictionsFromDescFile(const Graph &g,
                        const omnigraph::GraphElementFinder<Graph> &element_finder,
                        const std::string &desc_file) {
    fs::CheckFileExistenceFATAL(desc_file);
    std::ifstream descs(desc_file);
    std::string l;
    PartialGeneInfo info;

    size_t i = 1;
    GeneInitSeq starting_seqs;
    while (std::getline(descs, l)) {
        std::istringstream ss(l);
        ss >> info;

        INFO("Gene prediction #" << i << ": " << info);
        EdgeId e = element_finder.ReturnEdgeId(info.unitig_id);
        VERIFY(e != EdgeId());
        Range r = info.r;
        //exclusive right position
        r.end_pos += 1;
        //to 0-based coordinates
        r.shift(-1);
        VERIFY(r.size() > g.k());
        //to k+1-mer coordinates
//        r.end_pos -= g.k();

        if (!info.strand) {
            r = ConjugateRange(g, e, r);
            e = g.conjugate(e);
        }

        starting_seqs.insert(make_pair(info.gene_id, g.EdgeNucls(e).Subseq(r.start_pos, r.end_pos).str()));
        ++i;
    }
    return starting_seqs;
}

std::unordered_map<std::string, size_t>
CDSLengthsFromFile(const std::string &fn) {
    INFO("Parsing estimated CDS lengths from " << fn);
    fs::CheckFileExistenceFATAL(fn);
    std::unordered_map<std::string, size_t> answer;
    std::ifstream in(fn);
    std::string l;
    std::string gene_id;
    double est_len;
    while (std::getline(in, l)) {
        std::istringstream ss(l);
        ss >> gene_id;
        ss >> est_len;
        answer[gene_id] = RoundedProduct(3, est_len);
    }
    return answer;
}

static std::string GeneNameFromFasta(const std::string &header) {
    std::stringstream ss(header);
    std::string s;
    ss >> s;
    ss >> s;
    return s;
}

static GeneInitSeq PredictionsFromFastaFile(const std::string &fasta_fn) {
    fs::CheckFileExistenceFATAL(fasta_fn);
    io::FileReadStream gene_frags(fasta_fn);

    io::SingleRead gene_frag;

    size_t i = 1;
    GeneInitSeq starting_seqs;
    while (!gene_frags.eof()) {
        gene_frags >> gene_frag;
        INFO("Gene prediction #" << i << ": " << gene_frag.name());
        starting_seqs.insert(make_pair(GeneNameFromFasta(gene_frag.name()), gene_frag.GetSequenceString()));
        ++i;
    }
    return starting_seqs;
}

static void WriteComponent(const GraphComponent<Graph> &component, const std::string &prefix,
                           const std::set<GraphPos> &stop_codon_poss, const io::EdgeNamingF<Graph> &naming_f) {

    const auto &g = component.g();
    subgraph_extraction::WriteComponentWithDeadends(component, prefix, naming_f);

    INFO("Writing potential stop-codon positions to " << prefix << ".stops")
    std::ofstream stop_codon_os(prefix + ".stops");
    io::CanonicalEdgeHelper<Graph> canonical_helper(g, naming_f);
//1-based coordinate gives the start of the stop-codon
    for (GraphPos gpos : stop_codon_poss) {
        if (component.edges().count(gpos.first)) {
            stop_codon_os << canonical_helper.EdgeOrientationString(gpos.first, "\t")
                          << "\t" << (gpos.second + g.k() - 2) << "\n";
        } else {
            WARN("Earlier detected stop codon " << g.str(gpos.first) << " "
                                                << gpos.second << " is outside the component");
        }
    }
}

int main(int argc, char** argv) {
    utils::segfault_handler sh;
    gcfg cfg;

    process_cmdline(argc, argv, cfg);

    toolchain::create_console_logger();
    START_BANNER("SPAdes gene subgraph extractor");

    try {
        unsigned nthreads = cfg.nthreads;
        unsigned k = cfg.k;
        std::string out_folder = cfg.outdir + "/";
        fs::make_dirs(out_folder);

        INFO("K-mer length set to " << k);

        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);
        INFO("# of threads to use: " << nthreads);

        std::string tmpdir = cfg.tmpdir.empty() ? out_folder + "tmp" : cfg.tmpdir;
        fs::make_dirs(tmpdir);
        conj_graph_pack gp(k, tmpdir, 0);

        omnigraph::GraphElementFinder<Graph> element_finder(gp.g);
        INFO("Loading de Bruijn graph from " << cfg.graph);
        gp.kmer_mapper.Attach(); // TODO unnecessary
        io::EdgeLabelHelper<Graph> label_helper(element_finder,
                toolchain::LoadGraph(gp, cfg.graph));

        gp.EnsureBasicMapping();
        auto mapper = MapperInstance(gp);

        VERIFY(cfg.genes_desc.empty() != cfg.genes_seq.empty());

        auto starting_seqs = cfg.genes_desc.empty() ?
                              PredictionsFromFastaFile(cfg.genes_seq) :
                              PredictionsFromDescFile(gp.g, element_finder, cfg.genes_desc);

        const Graph &g = gp.g;
        PartialGeneProcessor processor(g, label_helper.edge_naming_f());
        auto cds_len_ests = CDSLengthsFromFile(cfg.cds_len_fn);

        const bool parallel = false;
        if (parallel) {
        INFO("Searching relevant subgraphs in parallel for all partial predictions");
        std::vector<std::string> flattened_ids;
        std::vector<std::string> flattened_part_genes;
        for (const auto &id_part : starting_seqs) {
            flattened_ids.push_back(id_part.first);
            flattened_part_genes.push_back(id_part.second);
        }

        size_t n = flattened_ids.size();

        std::vector<std::set<EdgeId>> flattened_relevant_edges(n);
        std::vector<std::set<GraphPos>> flattened_stop_poss(n);

        #pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < n; ++i) {
            flattened_relevant_edges[i] = processor.ProcessPartialGene(flattened_part_genes[i], mapper,
                                                                       utils::get(cds_len_ests, flattened_ids[i]),
                                                                       flattened_stop_poss[i]).edges();
        }
        INFO("Done searching subgraphs");

        std::set<EdgeId> edges;
        std::set<GraphPos> stop_codon_poss;
        for (size_t i = 0; i < n; ++i) {
            const std::string &gene_id = flattened_ids[i];
            utils::insert_all(edges, flattened_relevant_edges[i]);
            utils::insert_all(stop_codon_poss, flattened_stop_poss[i]);

            if (i == n - 1 || flattened_ids[i + 1] != gene_id) {
                subgraph_extraction::ComponentExpander expander(gp.g);
                auto component = expander.Expand(GraphComponent<Graph>::FromEdges(gp.g, edges.begin(),
                                                                                  edges.end(), /*add conjugate*/true));

                if (component.e_size() > 0) {
                    WriteComponent(component, out_folder + gene_id, stop_codon_poss, label_helper.edge_naming_f());
                } else {
                    INFO("Couldn't find a non-trivial component for gene " << gene_id);
                }
            }
        }
        } else {
        for (const auto &gene_id : utils::key_set(starting_seqs)) {
            INFO("Processing gene " << gene_id);
            INFO("Subgraphs extracted for all " << starting_seqs.count(gene_id) << " of its partial predictions will be united");
            std::set<EdgeId> edges;
            std::set<GraphPos> stop_codon_poss;
            for (const std::string &s : utils::get_all(starting_seqs, gene_id)) {
                utils::insert_all(edges, processor.ProcessPartialGene(s, mapper,
                                                                      utils::get(cds_len_ests, gene_id),
                                                                      stop_codon_poss).edges());
            }
            subgraph_extraction::ComponentExpander expander(gp.g);
            auto component = expander.Expand(GraphComponent<Graph>::FromEdges(gp.g, edges.begin(),
                                                                              edges.end(), /*add conjugate*/true));

            if (component.e_size() > 0) {
                WriteComponent(component, out_folder + gene_id, stop_codon_poss, label_helper.edge_naming_f());
            } else {
                INFO("Couldn't find a non-trivial component for gene " << gene_id);
            }
        }
        }

        INFO("Done");
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }
}
