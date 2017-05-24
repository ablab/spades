#include "modules/genome_consistance_checker.hpp"
#include "modules/path_extend/paired_library.hpp"
#include "assembly_graph/core/graph.hpp"
#include <algorithm>
#include <numeric>
#include <limits>

namespace debruijn_graph {
using omnigraph::MappingRange;
using namespace std;

//gap or overlap size. WITHOUT SIGN!
size_t AbsGap(const Range &a, const Range &b) {
    return max(a.end_pos, b.start_pos) - min(a.end_pos, b.start_pos);
}

bool GenomeConsistenceChecker::Consequent(const MappingRange &mr1, const MappingRange &mr2) const {
    //do not want to think about handling gaps near 0 position.
    size_t max_gap = max(AbsGap(mr1.initial_range, mr2.initial_range),
                         AbsGap(mr1.mapped_range, mr2.mapped_range));

    if (max_gap > absolute_max_gap_)
        return false;
    size_t len = max(min(mr1.initial_range.size(), mr1.mapped_range.size()),
                     min(mr2.initial_range.size(), mr2.mapped_range.size()));
    return max_gap <= size_t(math::round(relative_max_gap_* double(len)));
}

PathScore GenomeConsistenceChecker::CountMisassemblies(const BidirectionalPath &path) const {
    PathScore score = InternalCountMisassemblies(path);
    if (path.Size() == 0) {
        WARN ("0 length path in GCChecker!!!");
        return PathScore(0,0,0);
    }
    size_t total_length = path.LengthAt(0);
//TODO: constant;
    if (total_length > score.mapped_length * 2) {
        if (total_length > SIGNIFICANT_LENGTH_LOWER_LIMIT) {
            INFO ("For path length " << total_length <<" mapped less than half of the path, skipping");
        }
        return PathScore(0,0,0);
    } else {
        return score;
    }
}

MappingPath<EdgeId> GenomeConsistenceChecker::ConstructEdgeOrder(const string& chr_name) const {
    vector<pair<EdgeId, MappingRange>> to_sort;
    DEBUG ("constructing edge order for chr " << chr_name);
    for (auto e: storage_) {
        set<MappingRange> mappings = gp_.edge_pos.GetEdgePositions(e, chr_name);
        VERIFY_MSG(mappings.size() <= 1, "Presumably unique edge " << e << " with multiple mappings!");
        if (!mappings.empty()) {
            to_sort.push_back(make_pair(e, *mappings.begin()));
        }
    }
    DEBUG("Sorting " << to_sort << " positions:");
    sort(to_sort.begin(), to_sort.end(),
         [](const pair<EdgeId, MappingRange> & a, const pair<EdgeId, MappingRange> & b) {
        return a.second.initial_range.start_pos < b.second.initial_range.start_pos;
    });
    return MappingPathT(to_sort);
}

void GenomeConsistenceChecker::ReportEdge(EdgeId e, double w) const{
    INFO("Edge " << gp_.g.int_id(e) << " weight " << w << " len " << gp_.g.length(e) << " cov " << gp_.g.coverage(e));
    if (!genome_info_.Multiplicity(e)) {
        INFO(" no chromosome position");
    } else {
        auto info = genome_info_.UniqueChromosomeIdx(e);
        INFO ("Chromosome " << info.first << " index " << info.second);
    }
}

void GenomeConsistenceChecker::ReportVariants(vector<pair<double, EdgeId>> &sorted_w) const {
    sort(sorted_w.rbegin(), sorted_w.rend());
    size_t count = 0;
    double additional_weight = 0;
    size_t reporting = 4;
    for (const auto pair: sorted_w) {
        if (count == 0) {
            INFO("First candidate:");
        }
        if (count < reporting) {
            ReportEdge(pair.second, pair.first);
        } else {
            additional_weight += pair.first;
        }
        count++;
    }
    if (reporting < sorted_w.size()) {
        INFO("Additional weight " << additional_weight << " of " << sorted_w.size() - reporting <<
             " candidates");
    }
    if (sorted_w.size() == 0) {
        INFO("No uniqueness info");
    }
}

void GenomeConsistenceChecker::ReportPathEndByLongLib(const path_extend::BidirectionalPathSet &covering_paths,
                                                      EdgeId current_edge) const {
    vector<pair<double, EdgeId>> sorted_w;
    for (const auto & cov_path: covering_paths) {
        double w = cov_path->GetWeight();
        map<EdgeId, double> next_weigths;
        if (math::gr(w, 1.0)) {
            for (size_t p_ind = 0; p_ind < cov_path->Size(); p_ind++) {
                if (cov_path->At(p_ind) == current_edge) {
                    for (size_t p_ind2  = p_ind + 1; p_ind2 < cov_path->Size(); p_ind2++) {
                        if (gp_.g.length(cov_path->At(p_ind2)) >= storage_.min_length() ) {
                            next_weigths[cov_path->At(p_ind2)] += w;
                        }
                    }
                    break;
                }
            }
        }
        for (const auto &p: next_weigths) {
            sorted_w.push_back(make_pair(p.second, p.first));
        }
    }
    INFO("Looking on long reads, last long edge: ");
    ReportVariants(sorted_w);
}

void GenomeConsistenceChecker::ReportPathEndByPairedLib(const shared_ptr<path_extend::PairedInfoLibrary> paired_lib,
                                                        EdgeId current_edge) const {
    vector<pair<double, EdgeId>> sorted_w;
    set<EdgeId> result;
    paired_lib->FindJumpEdges(current_edge, result, std::numeric_limits<int>::min(), std::numeric_limits<int>::max(),
                              storage_.min_length());
    for (const auto e: result) {
        double w = paired_lib->CountPairedInfo(current_edge, e, std::numeric_limits<int>::min(),
                                               std::numeric_limits<int>::max());
        if (math::gr(w, 1.0))
            sorted_w.push_back(make_pair(w, e));
    }
    INFO("Looking on lib IS " << paired_lib->GetIS());
    ReportVariants(sorted_w);
}

void GenomeConsistenceChecker::CheckPathEnd(const BidirectionalPath &path) const {
    for (int i =  (int)path.Size() - 1; i >= 0; --i) {
        if (storage_.IsUnique(path.At(i))) {
            EdgeId current_edge = path.At(i);
            if (genome_info_.Multiplicity(current_edge)) {
                const auto &chr_info = genome_info_.UniqueChromosomeInfo(current_edge);
                size_t index = chr_info.UniqueEdgeIdx(current_edge);
                if (index == 0 || index == chr_info.size()) {
                    DEBUG("Path length " << path.Length() << " ended at the chromosome " << chr_info.name()
                          << (index == 0 ? " start": " end"));
                    return;
                }
            }
            INFO("Path length " << path.Length() << " ended, last unique: ");
            ReportEdge(current_edge, -1.0);
            for (size_t lib_index = 0; lib_index < reads_.lib_count(); ++lib_index) {
                const auto &lib = reads_[lib_index];
                if (lib.is_paired()) {
                    shared_ptr<path_extend::PairedInfoLibrary> paired_lib;
                    if (lib.is_mate_pair())
                        paired_lib = path_extend::MakeNewLib(gp_.g, lib, gp_.paired_indices[lib_index]);
                    else if (lib.type() == io::LibraryType::PairedEnd)
                        paired_lib = path_extend::MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
                    ReportPathEndByPairedLib(paired_lib, current_edge);
                } else if (lib.is_long_read_lib()) {
                    ReportPathEndByLongLib(long_reads_cov_map_[lib_index].GetCoveringPaths(current_edge), current_edge);
                }
            }
            return;
        }
    }
}

size_t GenomeConsistenceChecker::GetSupportingPathCount(EdgeId e1, EdgeId e2, size_t lib_index) const {
    auto covering_paths = long_reads_cov_map_[lib_index].GetCoveringPaths(e1);
    size_t res = 0;
    for (const auto & cov_path: covering_paths) {
        double w = cov_path->GetWeight();
        if (math::gr(w, 1.0)) {
            for (size_t p_ind = 0; p_ind < cov_path->Size(); p_ind++) {
                if (cov_path->At(p_ind) == e1) {
                    for (size_t p_ind2 = p_ind + 1; p_ind2 < cov_path->Size(); p_ind2++) {
                        if (storage_.IsUnique(cov_path->At(p_ind2))) {
                            if (e2 == cov_path->At(p_ind2))
                                res += size_t(w);
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }
    return res;
}

void GenomeConsistenceChecker::PrintMisassemblyInfo(EdgeId e1, EdgeId e2) const {
    VERIFY(genome_info_.Multiplicity(e1));
    VERIFY(genome_info_.Multiplicity(e2));
    const auto &chr_info1 = genome_info_.UniqueChromosomeInfo(e1);
    const auto &chr_info2 = genome_info_.UniqueChromosomeInfo(e2);
    size_t ind1 = chr_info1.UniqueEdgeIdx(e1);
    size_t ind2 = chr_info2.UniqueEdgeIdx(e2);
//FIXME: checks, compliment_strands;
    EdgeId true_next = chr_info1.EdgeAt((chr_info1.UniqueEdgeIdx(e1) + 1) % chr_info1.size());
    EdgeId true_prev = chr_info2.EdgeAt((chr_info2.UniqueEdgeIdx(e2) + chr_info2.size() - 1) % chr_info2.size());
    INFO("Next genomic edge " << true_next.int_id() << " len " << gp_.g.length(true_next) << " prev " << true_prev.int_id() << " len " << gp_.g.length(true_prev));
    if (chr_info1.name() == chr_info2.name() && ind1 < ind2) {
        INFO("Same chromosome large forward jump misassembly");
    } else if (chr_info1.name() == chr_info2.name() && ind1 > ind2)  {
        INFO("Backward jump misassembly");
    } else if (chr_info1.name().substr(1) == chr_info2.name().substr(1))  {
        string revers = (ind1 + ind2 + 2 > chr_info1.size() ? " backwards " : " forward " );
        INFO("Inversion" + revers +  "misassembly, chr edge size  " << chr_info1.size());
    } else if (ind1 + 1 == chr_info1.size()  || ind2 == 0) {
        string start_end = (ind2 == 0 ? " start " : " end ");
        INFO("Chromosome " + start_end + " misassembly ");
    } else {
        INFO("Something else misassembly");
    }
    for (size_t lib_index = 0; lib_index < reads_.lib_count(); ++lib_index) {
        const auto &lib = reads_[lib_index];
        if (lib.is_paired()) {
            shared_ptr<path_extend::PairedInfoLibrary> paired_lib;
            if (lib.is_mate_pair())
                paired_lib = path_extend::MakeNewLib(gp_.g, lib, gp_.paired_indices[lib_index]);
            else if (lib.type() == io::LibraryType::PairedEnd)
                paired_lib = path_extend::MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
            INFO("for lib " << lib_index << " IS" << paired_lib->GetIS());
            INFO("Misassembly weight regardless of dists: " << paired_lib->CountPairedInfo(e1, e2, -1000000, 1000000));
            INFO("Next weight " << paired_lib->CountPairedInfo(e1, true_next, -1000000, 1000000));
            INFO("Prev weight " << paired_lib->CountPairedInfo(true_prev, e2, -1000000, 1000000));
        } else if (lib.is_long_read_lib()) {
            INFO("for lib " << lib_index << " of long reads: ");
            INFO("Misassembly weight " << GetSupportingPathCount(e1, e2 ,lib_index));
            INFO("Next weight " << GetSupportingPathCount(e1, true_next ,lib_index) );
            INFO("Prev weight " << GetSupportingPathCount(true_prev, e2 ,lib_index) );

        }
    }
}

void GenomeConsistenceChecker::ClassifyPosition(size_t prev_pos, size_t cur_pos,
                                                const BidirectionalPath & path, PathScore &res) const{
    EdgeId cur_e = path.At(cur_pos);
    const auto& chr_info = genome_info_.UniqueChromosomeInfo(cur_e);
    size_t cur_in_genome = chr_info.UniqueEdgeIdx(cur_e);
    string cur_chr = chr_info.name();
    MappingRange cur_range = gp_.edge_pos.GetUniqueEdgePosition(cur_e, cur_chr);
    EdgeId prev_e = path.At(prev_pos);
    const auto& prev_chr_info = genome_info_.UniqueChromosomeInfo(prev_e);
    size_t prev_in_genome = prev_chr_info.UniqueEdgeIdx(prev_e);
    string prev_chr = prev_chr_info.name();
    MappingRange prev_range = gp_.edge_pos.GetUniqueEdgePosition(prev_e, prev_chr);

    res.mapped_length += cur_range.mapped_range.size();
    if (cur_in_genome == prev_in_genome + 1 && cur_chr == prev_chr) {
        int dist_in_genome = (int) cur_range.initial_range.start_pos -  (int) prev_range.initial_range.end_pos;
        int dist_in_path = (int) path.LengthAt(prev_pos) - (int) path.LengthAt(cur_pos) +
                (int) cur_range.mapped_range.start_pos - (int) prev_range.mapped_range.end_pos;
        DEBUG("Edge " << prev_e.int_id() << "  position in genome ordering: " << prev_in_genome);
        DEBUG("Gap in genome / gap in path: " << dist_in_genome << " / " << dist_in_path);
        if (size_t(abs(dist_in_genome - dist_in_path)) > absolute_max_gap_ &&
                (dist_in_genome * (1 + relative_max_gap_) < dist_in_path ||
                        dist_in_path * (1 + relative_max_gap_) < dist_in_genome)) {
            res.wrong_gap_size ++;
        }
    } else {
        if (cur_chr == prev_chr && (circular_edges_.find(prev_e) != circular_edges_.end() ||
                                    circular_edges_.find(cur_e) != circular_edges_.end())) {
            INFO("Skipping fake(circular) misassembly");
        } else if (cur_in_genome > prev_in_genome && cur_chr == prev_chr
                && prev_range.initial_range.end_pos + SIGNIFICANT_LENGTH_LOWER_LIMIT > cur_range.initial_range.start_pos) {
            INFO("Local misassembly between edges: "<<prev_e.int_id() << " and " << cur_e.int_id());
            size_t total = 0;
            for (auto j = prev_in_genome + 1; j < cur_in_genome; j ++) {
                total += gp_.g.length(chr_info.EdgeAt(j));
            }
            INFO("Jumped over " << cur_in_genome - prev_in_genome - 1 << " uniques of total length: " << total);
        } else if (IsCloseToEnd(prev_range, prev_chr_info) && IsCloseToStart(cur_range, chr_info)) {
            INFO ("Skipping fake misassembly - connected " << prev_chr << " and " << cur_chr);
        } else {
            INFO("Extensive misassembly between edges: "<<prev_e.int_id() << " and " << cur_e.int_id());
            INFO("Ranges: " << prev_range << " and " << cur_range);
            INFO("Genomic positions: " << prev_in_genome << ", " << prev_chr <<
                         " and " << cur_in_genome <<", "<< cur_chr<< " resp.");
            PrintMisassemblyInfo(prev_e, cur_e);
            res.misassemblies++;
        }
    }
}


PathScore GenomeConsistenceChecker::InternalCountMisassemblies(const BidirectionalPath &path) const {
    PathScore res(0, 0, 0);
    size_t prev_pos = std::numeric_limits<std::size_t>::max();
    for (int i = 0; i < (int) path.Size(); i++) {
//const method, so at instead of []
        EdgeId e = path.At(i);
        if (genome_info_.Multiplicity(e)) {
            if (prev_pos != std::numeric_limits<std::size_t>::max()) {
                ClassifyPosition(prev_pos, i, path, res);
            }
            prev_pos = i;
        }
    }
    return res;
}

vector<MappingRange> GenomeConsistenceChecker::FindBestRangeSequence(const set<MappingRange>& mappings) const {
    vector<MappingRange> to_process(mappings.begin(), mappings.end());
    sort(to_process.begin(), to_process.end(), [](const MappingRange & a, const MappingRange & b) -> bool
    {
        return a.mapped_range.start_pos < b.mapped_range.start_pos;
    } );
    size_t sz = to_process.size();
//max weight path in orgraph of mappings
    TRACE("constructing mapping graph" << sz << " vertices");
    vector<vector<size_t>> consecutive_mappings(sz);
    for (size_t i = 0; i < sz; i++) {
        for (size_t j = i + 1; j < sz; j++) {
            if (Consequent(to_process[i], to_process[j])) {
                consecutive_mappings[i].push_back(j);
            } else {
                if (to_process[j].mapped_range.start_pos > to_process[i].mapped_range.end_pos + absolute_max_gap_) {
                    break;
                }
            }
        }
    }
    vector<size_t> scores(sz), prev(sz);
    for (size_t i = 0; i < sz; i++) {
        scores[i] = to_process[i].initial_range.size();
        prev[i] = std::numeric_limits<std::size_t>::max();
    }
    for (size_t i = 0; i < sz; i++) {
        for (size_t j = 0; j < consecutive_mappings[i].size(); j++) {
            TRACE(consecutive_mappings[i][j]);
            if (scores[consecutive_mappings[i][j]] < scores[i]
                                                     + to_process[consecutive_mappings[i][j]].initial_range.size()) {
                scores[consecutive_mappings[i][j]] = scores[i]
                                                     + to_process[consecutive_mappings[i][j]].initial_range.size();
                prev[consecutive_mappings[i][j]] = i;
            }
        }
    }
    size_t cur_max = 0;
    size_t cur_i = 0;
    for (size_t i = 0; i < sz; i++) {
        if (scores[i] > cur_max) {
            cur_max = scores[i];
            cur_i = i;
        }
    }

    vector<MappingRange> answer;
    while (cur_i != std::numeric_limits<std::size_t>::max()) {
        answer.push_back(to_process[cur_i]);
        cur_i = prev[cur_i];
    }
    reverse(answer.begin(), answer.end());
    return answer;
}

map<EdgeId, string> GenomeConsistenceChecker::EdgeLabels() const {
    INFO("Constructing reference labels");
    map<EdgeId, string> answer;
    size_t count = 0;
    for (const auto &chr: genome_info_.Chromosomes()) {
        const auto &chr_info = genome_info_.ChrInfo(chr);
        for (size_t pos = 0; pos < chr_info.size(); ++pos) {
            EdgeId e = chr_info.EdgeAt(pos);
            auto mr = gp_.edge_pos.GetUniqueEdgePosition(e, chr);
            VERIFY(!answer.count(e));
            answer[e] += chr +
                         "order: " + to_string(count) +
                         "\n mapped range: " +
                         to_string(mr.mapped_range.start_pos) + " : "
                         + to_string(mr.mapped_range.end_pos) +
                         "\n init range: " +
                         to_string(mr.initial_range.start_pos) + " : "
                         + to_string(mr.initial_range.end_pos) + "\n";
        }
    }
    return answer;
}

void GenomeConsistenceChecker::Fill() {
    gp_.edge_pos.clear();
    if (!gp_.edge_pos.IsAttached()) {
        gp_.edge_pos.Attach();
    }

    //FIXME set the parameters to something more reasonable
    EdgesPositionHandler<Graph> tmp_edge_pos(gp_.g, 0, 0);
    visualization::position_filler::PosFiller<Graph> pos_filler(gp_.g, MapperInstance(gp_), tmp_edge_pos);

    for (const auto &chr: gp_.genome.GetChromosomes()) {
        pos_filler.Process(chr.sequence, "0_" + chr.name);
        pos_filler.Process(ReverseComplement(chr.sequence), "1_" + chr.name);
    }

    for (auto e: storage_) {
        FillPos(e, tmp_edge_pos);
    }

    vector<size_t> theoretic_lens;
    for (const auto &prefix: vector<std::string>{"0_", "1_"}) {
        for (const auto &chr: gp_.genome.GetChromosomes()) {
            string label = prefix + chr.name;
            INFO("Spelling label " << label);
            auto mapping_path = ConstructEdgeOrder(label);
            genome_info_.AddInfo(ChromosomeInfo(label, mapping_path));
            utils::push_back_all(theoretic_lens, MappedRegions(mapping_path));
        }
    }

    TheoreticLenStats(theoretic_lens);
}

void GenomeConsistenceChecker::TheoreticLenStats(vector<size_t> theoretic_lens) const {
    size_t total_len = std::accumulate(theoretic_lens.begin(), theoretic_lens.end(),
                                       0, std::plus<size_t>());

    std::sort(theoretic_lens.begin(), theoretic_lens.end());
    std::reverse(theoretic_lens.begin(), theoretic_lens.end());
    size_t cur = 0;
    size_t i = 0;
    while (cur < total_len / 2) {
        cur += theoretic_lens[i];
        i++;
    }
    INFO("Assuming gaps of length > " << storage_.min_length() << " unresolvable..");
    if (theoretic_lens.size() > 0)
        INFO("Rough estimates on N50/L50:" << theoretic_lens[i - 1] << " / " << i - 1 << " with len " << total_len);
}

map<string, size_t>
GenomeConsistenceChecker::TotalAlignedLengths(const EdgesPositionHandler<Graph> &tmp_edge_pos, EdgeId e) const {
    map<string, size_t> chr2len;
    for (const auto &edge_pos: tmp_edge_pos.GetEdgePositions(e)) {
        chr2len[edge_pos.contigId] += edge_pos.mr.initial_range.size();
    }
    return chr2len;
}

vector<size_t> GenomeConsistenceChecker::MappedRegions(const GenomeConsistenceChecker::MappingPathT &mapping_path) const {
    vector<size_t> mapped_regions;
    if (mapping_path.size() == 0)
        return mapped_regions;
    size_t pos = mapping_path.front().second.initial_range.start_pos;
    for (size_t i = 0; i < mapping_path.size(); i++) {
        auto current_range = mapping_path[i].second;
        INFO("Pos: " << i << " init_range " << current_range.initial_range
                     << " mapped to edge " << gp_.g.str(mapping_path[i].first)
                     << " range " << current_range.mapped_range);

        size_t curr_start = current_range.initial_range.start_pos;
        if (i > 0) {
            auto prev_range = mapping_path[i - 1].second;
            size_t prev_end = prev_range.initial_range.end_pos;
            if (curr_start - prev_end > unresolvable_len_) {
                INFO ("Large gap " << current_range.initial_range.start_pos -
                                      prev_range.initial_range.end_pos);
                mapped_regions.push_back(prev_end - pos);
                pos = curr_start;
            }
        }
    }
    mapped_regions.push_back(mapping_path.back().second.initial_range.end_pos - pos);
    return mapped_regions;
}

void GenomeConsistenceChecker::FillPos(EdgeId e, const EdgesPositionHandler<Graph> &tmp_edge_pos) {
    size_t total_mapped;
    string chr = ChromosomeByUniqueEdge(e, tmp_edge_pos, total_mapped);
    if (chr.empty())
        return;

    auto mapping_info = Merge(FindBestRangeSequence(tmp_edge_pos.GetEdgePositions(e, chr)));

    //FIXME what is the logic here?
    //used less that 0.9 of aligned length
    VERIFY(total_mapped >= mapping_info.second);
    if ((total_mapped - mapping_info.second) * 10 >=  gp_.g.length(e)) {
        INFO ("Edge " << gp_.g.int_id(e) << " length "<< gp_.g.length(e)  << "is potentially misassembled! mappings: ");
        for (auto mp : tmp_edge_pos.GetEdgePositions(e, chr)) {
            INFO("mp_range "<< mp.mapped_range.start_pos << " - " << mp.mapped_range.end_pos
                 << " init_range " << mp.initial_range.start_pos << " - " << mp.initial_range.end_pos );
            if (mp.initial_range.start_pos < absolute_max_gap_) {
                INFO ("Fake(linear order) misassembly on edge "<< e.int_id());
                circular_edges_.insert(e);
            }
        }
    }
    gp_.edge_pos.AddEdgePosition(e, chr, mapping_info.first);
}

pair<MappingRange, size_t> GenomeConsistenceChecker::Merge(const vector<MappingRange> &mappings) const {
    VERIFY(mappings.size() > 0);

    MappingRange mr = mappings.front();
    size_t total_mapped = mr.initial_range.size();
    for (size_t i = 1; i < mappings.size(); ++i) {
        total_mapped += mappings[i].initial_range.size();
        //FIXME why do we need merge?
        mr = mr.Merge(mappings[i]);
    }
    return make_pair(mr, total_mapped);
}

string GenomeConsistenceChecker::ChromosomeByUniqueEdge(const EdgeId &e,
                                                        const EdgesPositionHandler<Graph> &tmp_edge_pos,
                                                        size_t &total) const {
    DEBUG("Positioning edge " << gp_.g.str(e));
    map<string, size_t> total_al_lens = TotalAlignedLengths(tmp_edge_pos, e);
    total = 0;
    for (size_t c : utils::value_set(total_al_lens))
        total += c;

    if (total > size_t(math::round((double) gp_.g.length(e) * 1.5))) {
        INFO("Edge " << gp_.g.int_id(e) <<" was not unique due to the references, excluding ");
        return "";
    }

    string chr = "";
    size_t max_l = 0;
    for (const auto &p : total_al_lens) {
        if (p.second > max_l) {
            max_l = p.second;
            chr = p.first;
        }
    }

    DEBUG("Most likely chromosome " << chr << ". Mapped bp: " << max_l);
    //TODO: support non-unique edges;
    if (max_l < size_t(math::round((double) gp_.g.length(e) * 0.5))) {
        DEBUG("Too small a portion mapped. Edge not used");
        return "";
    }

    return chr;
};

}
