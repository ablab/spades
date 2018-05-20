#include "domain_graph_construction.hpp"

#include "assembly_graph/paths/path_processor.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"

#include "modules/path_extend/pe_utils.hpp"

#include "hmm/hmmfile.hpp"
#include "hmm/hmmmatcher.hpp"

#include "sequence/aa.hpp"
#include "io/reads/file_reader.hpp"
#include "io/reads/wrapper_collection.hpp"

#include "domain_graph.hpp"

namespace debruijn_graph {

//We should have file with .hmm or .hmm.gz extension
static std::string getFilename (const std::string &filename) {
    size_t found = filename.find_last_of("/");
    if (filename.substr(filename.size() - 3) == ".gz") {
        return filename.substr(found + 1, filename.size() - found - 8);
    } else if (filename.substr(filename.size() - 6) == ".fasta") {
        return filename.substr(found + 1, filename.size() - found - 7);
    } else {
        return filename.substr(found + 1, filename.size() - found - 5);
    }
}

template<class Graph>
class SetOfForbiddenEdgesPathChooser : public PathProcessor<Graph>::Callback {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph &g_;
    std::set<std::vector<EdgeId>> forbidden_edges_;
    std::vector<EdgeId> answer_path_;
    std::vector<size_t> peaks;

    bool CheckCoverageDiff(const path_extend::BidirectionalPath &path) const {
        double min_coverage = std::numeric_limits<double>::max();
        double max_coverage = std::numeric_limits<double>::min();

        for (size_t i = 0; i < path.Size(); ++i) {
            min_coverage = std::min(min_coverage, path.graph().coverage(path[i]));
            max_coverage = std::max(max_coverage, path.graph().coverage(path[i]));
        }
        return math::ge(50.0, max_coverage/min_coverage);
    }

    bool IsNewPathBetter(path_extend::BidirectionalPath &current, path_extend::BidirectionalPath &candidate) const {
        int current_length = (int)current.Length();
        int candidate_length = (int)candidate.Length();
        int diff_current = std::numeric_limits<int>::max();
        int diff_cand = std::numeric_limits<int>::max();
        for (auto peak : peaks) {
            diff_current = std::min(diff_current, abs((int)peak - current_length));
            diff_cand = std::min(diff_cand, abs((int)peak - candidate_length));
        }
        return diff_cand < diff_current;
    }

public:
    SetOfForbiddenEdgesPathChooser(const Graph &g, const std::set<std::vector<EdgeId>> &forbidden_edges) :
            g_(g), forbidden_edges_(forbidden_edges) {
    }

    void HandleReversedPath(const vector<EdgeId> &reversed_path) override {
        vector<EdgeId> forward_path = reversed_path;
        std::reverse(forward_path.begin(), forward_path.end());
        bool cross_start = false;
        bool cross_end = false;

        for (auto forbidden_path : forbidden_edges_) {
            if (forbidden_path.size() != 0) {
                if (find(std::begin(forward_path), std::end(forward_path), forbidden_path.front()) != forward_path.end()) {
                    cross_start = true;
                }
                if (find(std::begin(forward_path), std::end(forward_path), forbidden_path.back()) != forward_path.end()) {
                    cross_end = true;
                }
                if (cross_end && cross_start) {
                    return;
                }
            }
        }
        if (answer_path_.empty()) {
            if (!CheckCoverageDiff(path_extend::BidirectionalPath(g_, forward_path))) {
               return;
            }
            answer_path_ = forward_path;
            return;
        }

        if (!CheckCoverageDiff(path_extend::BidirectionalPath(g_,forward_path))) {
            return;
        }

        path_extend::BidirectionalPath current(g_, answer_path_);
        path_extend::BidirectionalPath candidate(g_, forward_path);
        if (IsNewPathBetter(current, candidate)) {
            answer_path_ = forward_path;
        }
    }

    void reset() {
        answer_path_.clear();
    }

    const vector<EdgeId>& answer() {
        return answer_path_;
    }

    void set_peaks(size_t first, size_t second) {
        peaks.clear();
        peaks.push_back(first);
        peaks.push_back(second);
    }
};

class DomainGraphConstructor {
public:
    DomainGraphConstructor(conj_graph_pack &gp)
            : gp_(gp) {}

    nrps::DomainGraph& ConstructGraph(const ContigAlnInfo &info) {
        gp_.EnsureIndex();
        INFO("A-domain graph construction started");
        ConstructNodes(info);
        ConstructStrongEdges();
        ConstructWeakEdges();
        INFO("A-domain graph construction ended");
        return graph;
    }

private:

    void ConnectWithWeakEdge(std::shared_ptr<nrps::Vertex> v1, std::shared_ptr<nrps::Vertex> v2,
                             SetOfForbiddenEdgesPathChooser<Graph> &chooser) {
        DEBUG("Trying to connect " << v1->name_ << " and " << v2->name_ << " with weak edge");
        int last_mapping = (int)gp_.g.length(mappings[v1->name_].back().first) - mappings[v1->name_].end_pos();
        int first_mapping = (int)mappings[v2->name_].start_pos();

        if (5000 < last_mapping + first_mapping) {
            return;
        }
        int min_len = 0;
        chooser.set_peaks(std::max(0, 500 - last_mapping - first_mapping), std::max(0, 1000 - last_mapping - first_mapping));
        if (gp_.g.EdgeEnd(v1->domain_edges_in_row_.back()) != gp_.g.EdgeStart(v2->domain_edges_in_row_[0])) {
            DEBUG("Trying to find paths from " << gp_.g.EdgeEnd(v1->domain_edges_in_row_.back()) << " to " << gp_.g.EdgeStart(v2->domain_edges_in_row_[0]));
            ProcessPaths(gp_.g, min_len, 4000 - last_mapping - first_mapping, gp_.g.EdgeEnd(v1->domain_edges_in_row_.back()), gp_.g.EdgeStart(v2->domain_edges_in_row_[0]), chooser);
            if(!chooser.answer().empty()) {
                DEBUG("Path was found");
                path_extend::BidirectionalPath p(gp_.g, chooser.answer());
                DEBUG("Start vertex: " << gp_.g.EdgeStart(p.Front()).int_id());
                DEBUG("End vertex: " << gp_.g.EdgeEnd(p.Back()).int_id());
                DEBUG("Path:");
                p.PrintDEBUG();
                graph.addEdge(v1, v2, false, p.Length() + last_mapping + first_mapping, chooser.answer());
            }
            else {
                DEBUG("Path was not found");
            }
        } else {
            graph.addEdge(v1, v2, false, last_mapping + first_mapping, vector<EdgeId>());
        }
        chooser.reset();
    }

    void ConstructWeakEdges() {
        std::set<std::vector<EdgeId>> forbidden_edges;
        for (const auto &mapping : mappings) {
            forbidden_edges.insert(mapping.second.simple_path());
        }
        SetOfForbiddenEdgesPathChooser<Graph> chooser(gp_.g, forbidden_edges);
        for (auto v1 : graph.getNodeSet()) {
            if (!graph.HasStrongEdge(v1) && (v1->near_to_the_end_of_contig_)) {
                for (auto v2 : graph.getNodeSet()) {
                    if (v1 != v2 && v1->rc_ != v2 && !graph.HasStrongIncomingEdge(v2) && v2->near_to_the_start_of_contig_) {
                        ConnectWithWeakEdge(v1, v2, chooser);
                    }
                }
            }
        }
    }

    class PairComparator
    {
        public:
            int operator()(const std::pair<int,int> &lhs, const std::pair<int,int> &rhs) const
            {
                return lhs.first < rhs.first;
            }
    };

    pair<int,int> SearchForSubvector(path_extend::BidirectionalPath *scaffold, MappingPath<EdgeId> &domain) {
        pair<int,int> answer;
        bool found = false;
        if (domain.size() > scaffold->Size()) {
            return std::make_pair<int,int>(-1,-1);
        }
        for (size_t i = 0; i < scaffold->Size() - domain.size() + 1; ++i) {
            if (found)
                break;
            for (size_t j = 0; j < domain.size(); ++j) {
                if ((*scaffold)[i + j] != domain[j].first) {
                    break;
                }
                if (j == domain.size() - 1) {
                    found = true;
                    answer.first = (int)i;
                    answer.second = (int)i+j;
                }
            }
        }
        if (!found) {
            return std::make_pair<int,int>(-1,-1);
        } else {
            return answer;
        }
    }

    pair<int,int> FindMappingToPath(path_extend::BidirectionalPath *scaffold, MappingPath<EdgeId> &domain, std::vector<EdgeId> &edges) {
        auto res = SearchForSubvector(scaffold, domain);
        if (res.first == -1) {
            return std::make_pair<int,int>(-1,-1);
        }
        int start = 0;
        size_t index = 0;
        for (;index != res.first; ++index) {
            start += gp_.g.length((*scaffold)[index]);
        }

        start += domain[0].second.mapped_range.start_pos;
        size_t end_index = index + domain.size() - 1;
        int end = 0;
        for (size_t i = 0; i < end_index; ++i) {
            end += gp_.g.length((*scaffold)[i]);
        }
        end += domain[domain.size() - 1].second.mapped_range.end_pos;
        for (int i = index; i <= end_index; ++i) {
            edges.push_back((*scaffold)[i]);
        }
        return std::make_pair(start, end);
    }

    void ConstructStrongEdgesInternal(std::pair<const std::string, MappingPath<EdgeId>> &mappings, path_extend::BidirectionalPath *path,
            std::map<size_t, std::map<std::pair<int, int>, std::pair<std::string, std::vector<EdgeId>>, PairComparator>> &mappings_for_path) {
        std::vector<EdgeId> edges;
        std::pair<int, int> coords = FindMappingToPath(path, mappings.second, edges);
        if (coords.first == -1) {
            return;
        }
        mappings_for_path[path->GetId()][coords] = std::make_pair(mappings.first, edges);

        if (coords.second + 5000 < path->Length()) {
           auto v = graph.getNode(mappings.first);
           v->near_to_the_end_of_contig_ = false;
        }
        if (coords.first > 5000) {
            auto v = graph.getNode(mappings.first);
            v->near_to_the_start_of_contig_ = false;
        }
    }

    size_t GetIndexFromPosition(size_t position, path_extend::BidirectionalPath *path) {
        size_t index = 0;
        for (index = 0; index < path->Size(); ++index) {
            size_t current_pos = path->LengthAt(index);
            if (current_pos > position) {
                break;
            }
        }
        return index;
    }

    vector<EdgeId> FindEdgesBetweenMappings(int first_mapping_end_coord, int second_mapping_start_coord, path_extend::BidirectionalPath *path) {
        vector<EdgeId> answer;
        if (first_mapping_end_coord < 0 || second_mapping_start_coord < 0) {
            return answer;
        }
        size_t first_mapping_end = GetIndexFromPosition(first_mapping_end_coord, path);
        size_t second_mapping_start = GetIndexFromPosition(second_mapping_start_coord, path);

        if (first_mapping_end > second_mapping_start) {
            return answer;
        }

        for (size_t i = first_mapping_end + 1; i < second_mapping_start; ++i) {
            answer.push_back((*path)[i]);
        }
        return answer;
    }

    void ConstructStrongEdges() {
        path_extend::GraphCoverageMap coverage_map(gp_.g, gp_.contig_paths);
        std::map<size_t, std::map<std::pair<int, int>, std::pair<std::string, std::vector<EdgeId>>, PairComparator>> mappings_for_path;
        std::map<size_t, path_extend::BidirectionalPath*> from_id_to_path;
        for (auto domain : mappings) {
            DEBUG("Processing mapping " << domain.first);
            auto mapping_path = domain.second;
            if (!mapping_path.empty()) {
                EdgeId first = mapping_path.front().first;
                auto path_container = coverage_map.GetCoveringPaths(first);
                for (auto path_pair : path_container) {
                    from_id_to_path[path_pair->GetId()] = path_pair;
                    from_id_to_path[path_pair->GetConjPath()->GetId()] = path_pair->GetConjPath();
                    ConstructStrongEdgesInternal(domain, path_pair, mappings_for_path);
                    ConstructStrongEdgesInternal(domain, path_pair->GetConjPath(), mappings_for_path);
                }
            }
        }
        for (auto p : mappings_for_path) {
            DEBUG("Processing path " << p.first);
            std::pair<std::pair<int, int>, std::pair<std::string, vector<EdgeId>>> prev(make_pair(-1, -1), make_pair("", std::vector<EdgeId>()));
            for (const auto& maps : p.second) {
                DEBUG("Processing mapping " << maps.second.first);
                DEBUG("Mapping start: " << maps.first.first << ". Mapping end: " << maps.first.second);

                if (prev.first.first == -1) {
                    prev = maps;
                    continue;
                }
                if (prev.first.second > maps.first.first) {
                    DEBUG("Mapping intersects with other, skipping");
                    graph.removeVertex(maps.second.first);
                    continue;
                }

                if (prev.first.second < maps.first.first && maps.first.first - prev.first.second < 20000) {
                    DEBUG("Connecting " << prev.second << " and " << maps.second);
                    graph.addEdge(prev.second.first, maps.second.first, true, maps.first.first - prev.first.second, FindEdgesBetweenMappings(prev.first.second, maps.first.first, from_id_to_path[p.first]));
                }
                prev = maps;
            }
        }
    }

    //TODO: try some good coverage strategy
    bool IsInsideRepeat(std::shared_ptr<nrps::Vertex> v) {
        if (v->domain_edges_in_row_.size() > 1) {
            return false;
        }
        EdgeId e = v->domain_edges_in_row_[0];
        if (gp_.g.IncomingEdgeCount(gp_.g.EdgeStart(e)) > 1 ||  gp_.g.OutgoingEdgeCount(gp_.g.EdgeEnd(e)) > 1) {
            return true;
        }
        return false;
    }

    void ConstructNodes(const ContigAlnInfo &info) {
        auto mapper = MapperInstance(gp_);

        unsigned id = 1;
        for (const auto &aln : info) {
            Sequence sequence = Sequence(aln.seq);
            std::string name = aln.name + "_" + std::to_string(id), name_rc = name + "rc";

            auto edges = mapper->MapSequence(sequence);
            if (edges.simple_path().size() == 0)
                continue;

            DEBUG("Adding vertex " << name);
            graph.addVertex(name,
                            edges.simple_path(),
                            edges.front().second.mapped_range.start_pos, edges.back().second.mapped_range.end_pos,
                            aln.type);
            auto rc_edges = mapper->MapSequence(!sequence);
            DEBUG("Adding vertex " << name_rc);
            graph.addVertex(name_rc,
                            rc_edges.simple_path(),
                            rc_edges.front().second.mapped_range.start_pos, rc_edges.back().second.mapped_range.end_pos,
                            aln.type);
            mappings[name] = edges;
            mappings[name_rc] = rc_edges;
            graph.makeRC(name, name_rc);
            id++;
        }

        for (auto v : graph.getVertexSet()) {
            v->max_visited_ = IsInsideRepeat(v) ? 2 : 1;
        }
    }

    conj_graph_pack &gp_;
    nrps::DomainGraph graph;
    std::map<std::string, MappingPath<EdgeId>> mappings;
    DECL_LOGGER("AGraph");
};

static std::vector<std::string> getFileVector(const std::string &hmm_files) {
    std::string s = hmm_files;
    std::vector<std::string> result;
    std::string delimiter = ",";
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        result.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    result.push_back(s);
    return result;
}

static void match_contigs_internal(hmmer::HMMMatcher &matcher, path_extend::BidirectionalPath* path,
                                   const std::string &path_string,
                                   const std::string &type, const hmmer::HMM &hmm, ContigAlnInfo &res, io::OFastaReadStream &oss_contig) {
    for (size_t shift = 0; shift < 3; ++shift) {
        std::string ref_shift = std::to_string(path->GetId()) + "_" + std::to_string(shift);
        std::string seq_aas = aa::translate(path_string.c_str() + shift);
        matcher.match(ref_shift.c_str(), seq_aas.c_str());
    }
    matcher.summarize();

    for (const auto &hit : matcher.hits()) {
        if (!hit.reported() || !hit.included())
            continue;

        for (const auto &domain : hit.domains()) {
            std::pair<int, int> seqpos = domain.seqpos();
            std::pair<int, int> seqpos2 = domain.hmmpos();
            INFO("First - " << seqpos2.first << ", second - " << seqpos2.second);
            INFO("First - " << seqpos.first << ", second - " << seqpos.second);

            if (seqpos.second - seqpos.first < hmm.length() / 2)
                continue;

            int shift = hit.name()[strlen(hit.name()) - 1] - '0';
            seqpos.first = seqpos.first * 3  + shift;
            seqpos.second = seqpos.second * 3  + shift;
            std::string name(hit.name());
            oss_contig << io::SingleRead(name, path_string);
            INFO(name);
            INFO("First - " << seqpos.first << ", second - " << seqpos.second);
            res.push_back({type, name, unsigned(seqpos.first), unsigned(seqpos.second), path_string.substr(seqpos.first, seqpos.second - seqpos.first)});
        }
    }
}

void match_contigs(const path_extend::PathContainer &contig_paths, const path_extend::ScaffoldSequenceMaker &scaffold_maker,
                   const std::string &type, const hmmer::HMM &hmm, const hmmer::hmmer_cfg &cfg,
                   ContigAlnInfo &res, io::OFastaReadStream &oss_contig) {
    INFO("Total contigs: " << contig_paths.size());
    INFO("Model length - " << hmm.length());
    for (auto iter = contig_paths.begin(); iter != contig_paths.end(); ++iter) {
        hmmer::HMMMatcher matcher(hmm, cfg);
        path_extend::BidirectionalPath* path = iter.get();
        if (path->Length() <= 0)
            continue;
        std::string path_string = scaffold_maker.MakeSequence(*path);
        match_contigs_internal(matcher, path, path_string, type, hmm, res, oss_contig);

        path_extend::BidirectionalPath* conj_path = path->GetConjPath();
        if (conj_path->Length() <= 0)
            continue;
        std::string path_string_conj = scaffold_maker.MakeSequence(*conj_path);
        match_contigs_internal(matcher, conj_path, path_string_conj, type, hmm, res, oss_contig);
    }
}

ContigAlnInfo DomainMatcher::MatchDomains(conj_graph_pack &gp) {
    if (fs::check_existence(cfg::get().output_dir + "/temp_anti/"))
        fs::remove_dir(cfg::get().output_dir + "/temp_anti/");
    fs::make_dirs(cfg::get().output_dir + "/temp_anti/");
    fs::make_dirs(cfg::get().output_dir + "/bgc_in_gfa/");

    ContigAlnInfo res;
    hmmer::hmmer_cfg hcfg;
    auto file_vector = getFileVector(cfg::get().hmm_set);
    path_extend::ScaffoldSequenceMaker scaffold_maker(gp.g);

    io::OFastaReadStream oss_contig(cfg::get().output_dir + "/temp_anti/restricted_edges.fasta");
    for (const auto &file : file_vector) {
        hmmer::HMMFile hmmfile(file);
        if (!hmmfile.valid())
        FATAL_ERROR("Error opening HMM file "<< file);
        auto hmmw = hmmfile.read();
        INFO("Matching contigs with " << file);
        match_contigs(gp.contig_paths, scaffold_maker,
                      getFilename(file), hmmw.get(), hcfg,
                      res, oss_contig);
    }

    return res;
}

void DomainGraphConstruction::run(conj_graph_pack &gp, const char*) {
    std::vector<std::string> domain_filenames;

    auto res = DomainMatcher().MatchDomains(gp);

    DomainGraphConstructor constructor(gp);
    auto &graph = constructor.ConstructGraph(res);

    graph.FindDomainOrderings(gp, cfg::get().output_dir + "/gene_clusters.fasta");
    graph.ExportToDot(cfg::get().output_dir + "/domain_graph.dot");
    INFO("Go to export paths");
    graph.ExportPaths(gp, cfg::get().output_dir);
}

}
