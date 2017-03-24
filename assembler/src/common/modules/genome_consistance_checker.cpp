#include "modules/genome_consistance_checker.hpp"
#include "modules/path_extend/paired_library.hpp"
#include "assembly_graph/core/graph.hpp"
#include <algorithm>
#include <limits>

namespace debruijn_graph {
using omnigraph::MappingRange;
using namespace std;

//gap or overlap size. WITHOUT SIGN!
static size_t gap(const Range &a, const Range &b) {
    return max(a.end_pos, b.start_pos) - min (a.end_pos, b.start_pos);
}
bool GenomeConsistenceChecker::consequent(const Range &mr1, const Range &mr2) const{
    if (mr1.end_pos > mr2.start_pos + absolute_max_gap_)
        return false;
    if (mr1.end_pos + absolute_max_gap_ < mr2.start_pos)
        return false;
    return true;

}
bool GenomeConsistenceChecker::consequent(const MappingRange &mr1, const MappingRange &mr2) const {
    //do not want to think about handling gaps near 0 position.
    if (!consequent(mr1.initial_range, mr2.initial_range) || !consequent(mr1.mapped_range, mr2.mapped_range))
        return false;
    size_t initial_gap = gap(mr1.initial_range, mr2.initial_range);
    size_t mapped_gap = gap(mr1.mapped_range, mr2.mapped_range);
    size_t max_gap = max(initial_gap, mapped_gap);
    if ( max_gap > relative_max_gap_* double (max (min(mr1.initial_range.size(), mr1.mapped_range.size()), min(mr2.initial_range.size(), mr2.mapped_range.size()))))
        return false;
    return true;
}

PathScore GenomeConsistenceChecker::CountMisassemblies(const BidirectionalPath &path) const {
    PathScore score = InternalCountMisassemblies(path);

    size_t total_length = path.LengthAt(0);
//TODO: constant;
    if (total_length > score.mapped_length * 2) {
        if (total_length > 10000) {
            INFO ("For path length " << total_length <<" mapped less than half of the path, skipping");
        }
        return PathScore(0,0,0);
    } else {
        return score;
    }
}

map<std::string, vector<pair<EdgeId, MappingRange>>> GenomeConsistenceChecker::ConstructEdgeOrder() const {
    auto chromosomes = gp_.genome.GetChromosomes();
    map<std::string, vector<pair<EdgeId, MappingRange>>> res;
    vector<std::string> fixed_prefixes{"0_", "1_"};
    for (const auto prefix: fixed_prefixes) {
        for (const auto &chr: chromosomes) {
            string label = prefix + chr.name;
            res[label] = ConstructEdgeOrder(label);
        }
    }
    return res;
}

vector<pair<EdgeId, MappingRange>> GenomeConsistenceChecker::ConstructEdgeOrder(const string chr_name) const {
    vector<pair<EdgeId, MappingRange> > to_sort;
    DEBUG ("constructing edge order for chr " << chr_name);
    for(auto e: storage_) {
        if (excluded_unique_.find(e) == excluded_unique_.end()) {
            set<MappingRange> mappings = gp_.edge_pos.GetEdgePositions(e, chr_name);
            if (mappings.size() > 1) {
                INFO("Presumably unique edge " << e << " with multiple mappings!");
            } else if (mappings.size() == 0) {
                continue;
            } else {
                to_sort.push_back(make_pair(e, *mappings.begin()));
            }
        }
    }
    DEBUG("Sorting " << to_sort << " positions:");
    sort(to_sort.begin(), to_sort.end(), [](const pair<EdgeId, MappingRange> & a, const pair<EdgeId, MappingRange> & b) -> bool
    {
        return a.second.initial_range.start_pos < b.second.initial_range.start_pos;
    }
    );
    return to_sort;
}

void GenomeConsistenceChecker::SpellChromosome(const pair<string, vector<pair<EdgeId, MappingRange>>>& chr_info, vector<size_t>& lengths) {
    INFO("Spelling label " << chr_info.first);
    chromosomes_spelled_[chr_info.first] = vector<EdgeId>(0);
    vector<size_t> starts;
    vector<size_t> ends;
    size_t count = 0;
    for(size_t i = 0; i < chr_info.second.size(); i++) {
        if (i > 0 && chr_info.second[i].second.initial_range.start_pos - chr_info.second[i-1].second.initial_range.end_pos > unresolvable_len_ ) {
            INFO ("Large gap " << chr_info.second[i].second.initial_range.start_pos - chr_info.second[i-1].second.initial_range.end_pos );
            starts.push_back(chr_info.second[i].second.initial_range.start_pos);
            ends.push_back(chr_info.second[i-1].second.initial_range.end_pos);
        }
        if (i == 0) {
            starts.push_back(chr_info.second[i].second.initial_range.start_pos);
        }
        if (i == chr_info.second.size() - 1){
            ends.push_back(chr_info.second[i].second.initial_range.end_pos);
        }
        INFO("Pos:" << count << " edge " << gp_.g.int_id(chr_info.second[i].first) << " length "<< gp_.g.length(chr_info.second[i].first) <<
             " coverage " << gp_.g.coverage(chr_info.second[i].first) << " mapped to " << chr_info.second[i].second.mapped_range.start_pos
             << " - " << chr_info.second[i].second.mapped_range.end_pos << " init_range " << chr_info.second[i].second.initial_range.start_pos << " - " << chr_info.second[i].second.initial_range.end_pos );
        genome_spelled_[chr_info.second[i].first] = std::make_pair(chr_info.first, count);
        chromosomes_spelled_[chr_info.first].push_back(chr_info.second[i].first);
        count++;
    }
    for (size_t i = 0; i < starts.size(); i++)
        lengths.push_back(ends[i] - starts[i]);

}

void GenomeConsistenceChecker::SpellGenome() {
    auto all_positions = ConstructEdgeOrder();
    vector<size_t> theoretical_scaff_lengths;
    for(const auto &chr_info: all_positions) {
        SpellChromosome(chr_info, theoretical_scaff_lengths);
    }
    size_t total_len = 0;
    for (size_t i = 0; i < theoretical_scaff_lengths.size(); i++) {
        total_len += theoretical_scaff_lengths[i];
    }
    sort(theoretical_scaff_lengths.begin(), theoretical_scaff_lengths.end());
    reverse(theoretical_scaff_lengths.begin(), theoretical_scaff_lengths.end());
    size_t cur = 0;
    size_t i = 0;
    while (cur < total_len / 2 && i < theoretical_scaff_lengths.size()) {
        cur += theoretical_scaff_lengths[i];
        i++;
    }
    INFO("Assuming gaps of length > " << storage_.GetMinLength() << " unresolvable..");
    if (theoretical_scaff_lengths.size() > 0)
        INFO("Rough estimates on N50/L50:" << theoretical_scaff_lengths[i - 1] << " / " << i - 1 << " with len " << total_len);
}

void GenomeConsistenceChecker::ReportEdge(EdgeId e, double w) const{
    INFO( "Edge " << gp_.g.int_id(e) << " weight " << w << " len " << gp_.g.length(e));
    if (genome_spelled_.find(e) != genome_spelled_.end()) {
        INFO ("Chromosome " << genome_spelled_.at(e).first << " index " << genome_spelled_.at(e).second );
    } else {
        INFO(" no chromosome position");
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
void GenomeConsistenceChecker::CheckPathEnd(const BidirectionalPath &path) const {
    for (int i =  (int)path.Size() - 1; i >= 0; --i) {
        if (storage_.IsUnique(path.At(i))) {
            EdgeId current_edge = path.At(i);
            auto dataset_info_ = cfg::get().ds;
            for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); ++lib_index) {
                const auto &lib = dataset_info_.reads[lib_index];
                vector<pair<double, EdgeId>> sorted_w;
                if (lib.is_paired()) {
                    shared_ptr<path_extend::PairedInfoLibrary> paired_lib;
                    if (lib.is_mate_pair())
                        paired_lib = path_extend::MakeNewLib(gp_.g, lib, gp_.paired_indices[lib_index]);
                    else if (lib.type() == io::LibraryType::PairedEnd)
                        paired_lib = path_extend::MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
                    set<EdgeId> result;
                    paired_lib->FindJumpEdges(current_edge, result, -1000000, 1000000, storage_.GetMinLength());
                    for (const auto e: result) {
                        double w = paired_lib->CountPairedInfo(current_edge, e, -1000000, 1000000);
                        if (math::gr(w, 1.0))
                            sorted_w.push_back(make_pair(w, e));
                    }
                    INFO("Path length " << path.Length() << "ended, looking on lib IS " <<
                         paired_lib->GetIS() << " last long edge: ");
                    ReportEdge(current_edge, 0.0);
                } else if (lib.is_long_read_lib()) {
                    auto covering_paths = long_reads_cov_map_[lib_index]->GetCoveringPaths(current_edge);
                    for (const auto & cov_path: covering_paths) {
                        double w = cov_path->GetWeight();
                        map<EdgeId, double> next_weigths;
                        if (math::gr(w, 1.0)) {
                            for (size_t p_ind = 0; p_ind < cov_path->Size(); p_ind++) {
                                if (cov_path->At(p_ind) == current_edge) {
                                    for (size_t p_ind2  = p_ind + 1; p_ind2 < cov_path->Size(); p_ind2++) {
                                        if (gp_.g.length(cov_path->At(p_ind2)) >= storage_.GetMinLength() ) {
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
                    INFO("Path length " << path.Length() << " looking on long reads lib " << lib_index << " last long edge: ");
                    ReportVariants(sorted_w);
                }

            }
            return;
        }
    }
}

size_t GenomeConsistenceChecker::GetSupportingPathCount(EdgeId e1, EdgeId e2, size_t lib_index) const {
    auto covering_paths = long_reads_cov_map_[lib_index]->GetCoveringPaths(e1);
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
    VERIFY(genome_spelled_.find(e1) != genome_spelled_.end());
    VERIFY(genome_spelled_.find(e2) != genome_spelled_.end());
    auto pair1 = genome_spelled_.at(e1);
    auto pair2 = genome_spelled_.at(e2);
//FIXME: checks, compliment_strands;
    size_t s1 = chromosomes_spelled_.at(pair1.first).size();
    size_t s2 = chromosomes_spelled_.at(pair2.first).size();
    EdgeId true_next = chromosomes_spelled_.at(pair1.first)[(pair1.second + 1) % s1];
    EdgeId true_prev = chromosomes_spelled_.at(pair2.first)[(pair2.second + s2 - 1) % s2];
    auto dataset_info_ = cfg::get().ds;
    for (size_t lib_index = 0; lib_index < dataset_info_.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info_.reads[lib_index];
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

PathScore GenomeConsistenceChecker::InternalCountMisassemblies(const BidirectionalPath &path) const {
    PathScore res(0, 0, 0);
    EdgeId prev;
    size_t prev_in_genome = std::numeric_limits<std::size_t>::max();
    size_t prev_in_path = std::numeric_limits<std::size_t>::max();
    string path_chr = "";
    MappingRange prev_range;
    for (int i = 0; i < (int) path.Size(); i++) {
        if (genome_spelled_.find(path.At(i)) != genome_spelled_.end()) {
//const method, so at instead of []
            size_t cur_in_genome =  genome_spelled_.at(path.At(i)).second;
            string cur_chr = genome_spelled_.at(path.At(i)).first;
            MappingRange cur_range = *gp_.edge_pos.GetEdgePositions(path.At(i), cur_chr).begin();
            if (prev_in_genome != std::numeric_limits<std::size_t>::max()) {
                if (cur_in_genome == prev_in_genome + 1 && cur_chr == path_chr) {
                    int dist_in_genome = (int) cur_range.initial_range.start_pos -  (int) prev_range.initial_range.end_pos;
                    int dist_in_path = (int) path.LengthAt(prev_in_path) - (int) path.LengthAt(i) +  (int) cur_range.mapped_range.start_pos - (int) prev_range.mapped_range.end_pos;
                    DEBUG("Edge " << prev.int_id() << "  position in genome ordering: " << prev_in_genome);
                    DEBUG("Gap in genome / gap in path: " << dist_in_genome << " / " << dist_in_path);
                    if (size_t(abs(dist_in_genome - dist_in_path)) > absolute_max_gap_ && (dist_in_genome * (1 + relative_max_gap_) < dist_in_path || dist_in_path * (1 + relative_max_gap_) < dist_in_genome)) {

                        res.wrong_gap_size ++;
                    }
                } else {
                    if (cur_chr != path_chr) {
                        DEBUG("interchromosome misassembly!");
                    }
                    if (cur_chr == path_chr && (circular_edges_.find(prev) != circular_edges_.end() ||
                                circular_edges_.find(path.At(i)) != circular_edges_.end())) {
                        INFO("Skipping fake(circular) misassembly");
                    } else {
                        INFO("Misassembly between edges: "<<prev.int_id() << " and " << path.At(i).int_id());
                        INFO("Ranges: " << prev_range << " and " << cur_range);
                        INFO("Genomic positions: " << prev_in_genome<< ", " << path_chr << " and " << cur_in_genome <<", "<< cur_chr<< " resp.");
                        PrintMisassemblyInfo(prev, path.At(i));
                        res.misassemblies++;
                    }
                }
            }
            res.mapped_length += cur_range.mapped_range.size();
            prev = path.At(i);
            prev_in_genome = cur_in_genome;
            prev_range = cur_range;
            path_chr = cur_chr;
            prev_in_path = i;
        }
    }
    if (prev_in_path != std::numeric_limits<std::size_t>::max())
        DEBUG("Edge " << prev.int_id() << "  position in genome ordering: " << prev_in_genome);
    return res;
}
void GenomeConsistenceChecker::RefillPos() {
    RefillPos("0");
    RefillPos("1");
}


void GenomeConsistenceChecker::RefillPos(const string &strand) {
    for (auto chr: gp_.genome.GetChromosomes()) {
        for (auto e: storage_) {
            RefillPos("_" + strand + "_" + chr.name, e);
        }
    }
}

void GenomeConsistenceChecker::FindBestRangeSequence(const set<MappingRange>& old_mappings, vector<MappingRange>& used_mappings) const {
    vector<MappingRange> to_process (old_mappings.begin(), old_mappings.end());
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
            if (consequent(to_process[i], to_process[j])) {
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
            if (scores[consecutive_mappings[i][j]] < scores[i] + to_process[consecutive_mappings[i][j]].initial_range.size()) {
                scores[consecutive_mappings[i][j]] = scores[i] + to_process[consecutive_mappings[i][j]].initial_range.size();
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
    used_mappings.clear();
    while (cur_i != std::numeric_limits<std::size_t>::max()) {
        used_mappings.push_back(to_process[cur_i]);
        cur_i = prev[cur_i];
    }
    reverse(used_mappings.begin(), used_mappings.end());
};

void GenomeConsistenceChecker::RefillPos(const string &strand, const EdgeId &e) {
    set<MappingRange> old_mappings = gp_.edge_pos.GetEdgePositions(e, strand);
    TRACE("old mappings sz " << old_mappings.size() );
    size_t total_mapped = 0;
    for (auto mp:old_mappings) {
        total_mapped += mp.initial_range.size();
    }
    if (total_mapped  > (double) gp_.g.length(e) * 1.5) {
       INFO ("Edge " << gp_.g.int_id(e) << "is not unique, excluding");
       excluded_unique_.insert(e);
       return;
    }
//TODO: support non-unique edges;
    if (total_mapped  < (double) gp_.g.length(e) * 0.5) {
        DEBUG ("Edge " << gp_.g.int_id(e) << "is not mapped on strand "<< strand <<", not used");
        return;
    }
    TRACE(total_mapped << " " << gp_.g.length(e));
    string new_strand = strand.substr(1);
    vector<MappingRange> used_mappings;
    FindBestRangeSequence(old_mappings, used_mappings);

    size_t cur_i = 0;
    MappingRange new_mapping;
    new_mapping = used_mappings[cur_i];
    size_t used_mapped = new_mapping.initial_range.size();
    TRACE ("Edge " << gp_.g.int_id(e) << " length "<< gp_.g.length(e));
    TRACE ("new_mapping mp_range "<< new_mapping.mapped_range.start_pos << " - " << new_mapping.mapped_range.end_pos
         << " init_range " << new_mapping.initial_range.start_pos << " - " << new_mapping.initial_range.end_pos );
    while (cur_i  < used_mappings.size() - 1) {
        cur_i ++;
        used_mapped += used_mappings[cur_i].initial_range.size();
        new_mapping = new_mapping.Merge(used_mappings[cur_i]);
        TRACE("new_mapping mp_range "<< new_mapping.mapped_range.start_pos << " - " << new_mapping.mapped_range.end_pos
             << " init_range " << new_mapping.initial_range.start_pos << " - " << new_mapping.initial_range.end_pos );
    }
//used less that 0.9 of aligned length
    if (total_mapped * 10  >=  used_mapped * 10  + gp_.g.length(e)) {
        INFO ("Edge " << gp_.g.int_id(e) << " length "<< gp_.g.length(e)  << "is potentially misassembled! mappings: ");
        for (auto mp:old_mappings) {
            INFO("mp_range "<< mp.mapped_range.start_pos << " - " << mp.mapped_range.end_pos << " init_range " << mp.initial_range.start_pos << " - " << mp.initial_range.end_pos );
            if (mp.initial_range.start_pos < absolute_max_gap_) {
                INFO ("Fake(linear order) misassembly on edge "<< e.int_id());
                circular_edges_.insert(e);
            }
        }

    }
    gp_.edge_pos.AddEdgePosition(e, new_strand, new_mapping);
}



}
