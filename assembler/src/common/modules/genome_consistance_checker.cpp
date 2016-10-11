#include "modules/genome_consistance_checker.hpp"
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
    PathScore straight = CountMisassembliesWithStrand(path, "0");
    PathScore reverse = CountMisassembliesWithStrand(path, "1");
    size_t total_length = path.LengthAt(0);
//TODO: constant;
    if (total_length > std::max(straight.mapped_length, reverse.mapped_length) * 2) {
        if (total_length > 10000) {
            INFO ("For path length " << total_length <<" mapped less than half of the path, skipping");
        }
        return PathScore(0,0,0);
    } else {
        if (straight.mapped_length > reverse.mapped_length) {
            return straight;
        } else {
            return reverse;
        }
    }
}

vector<pair<EdgeId, MappingRange> > GenomeConsistenceChecker::ConstructEdgeOrder() const {
    vector<pair<EdgeId, MappingRange> > to_sort;
    for(auto e: storage_) {
        if (excluded_unique_.find(e) == excluded_unique_.end() ) {
            set<MappingRange> mappings = gp_.edge_pos.GetEdgePositions(e, "fxd0");
            if (mappings.size() > 1) {
                INFO("edge " << e << "smth strange");
            } else if (mappings.size() == 0) {
                continue;
            } else {
                to_sort.push_back(make_pair(e, *mappings.begin()));
            }
        }
    }
    sort(to_sort.begin(), to_sort.end(), [](const pair<EdgeId, MappingRange> & a, const pair<EdgeId, MappingRange> & b) -> bool
    {
        return a.second.initial_range.start_pos < b.second.initial_range.start_pos;
    }
    );
    return to_sort;
}


void GenomeConsistenceChecker::SpellGenome() {
    size_t count = 0;
    auto to_sort = ConstructEdgeOrder();
    vector<size_t> starts;
    vector<size_t> ends;
    for(size_t i = 0; i <to_sort.size(); i++) {
        if (i > 0 && to_sort[i].second.initial_range.start_pos - to_sort[i-1].second.initial_range.end_pos > storage_.GetMinLength() ) {
            INFO ("Large gap " << to_sort[i].second.initial_range.start_pos - to_sort[i-1].second.initial_range.end_pos );
            starts.push_back(to_sort[i].second.initial_range.start_pos);
            ends.push_back(to_sort[i-1].second.initial_range.end_pos);
        }
        if (i == 0) {
            starts.push_back(to_sort[i].second.initial_range.start_pos);
        }
        if (i == to_sort.size() - 1){
            ends.push_back(to_sort[i].second.initial_range.end_pos);
        }
        INFO("edge " << gp_.g.int_id(to_sort[i].first) << " length "<< gp_.g.length(to_sort[i].first) <<
                     " coverage " << gp_.g.coverage(to_sort[i].first) << " mapped to " << to_sort[i].second.mapped_range.start_pos
             << " - " << to_sort[i].second.mapped_range.end_pos << " init_range " << to_sort[i].second.initial_range.start_pos << " - " << to_sort[i].second.initial_range.end_pos );
        genome_spelled_[to_sort[i].first] = count;
        count++;
    }
    vector<size_t> lengths;
    size_t total_len = 0;
    for (size_t i = 0; i < starts.size(); i++) {
        lengths.push_back(ends[i] - starts[i]);
        total_len += lengths[i];
    }
    sort(lengths.begin(), lengths.end());
    reverse(lengths.begin(), lengths.end());
    size_t cur = 0;
    size_t i = 0;
    while (cur < total_len / 2 && i < lengths.size()) {
        cur += lengths[i];
        i++;
    }
    INFO("Assuming gaps of length > " << storage_.GetMinLength() << " unresolvable..");
    if (lengths.size() > 0)
        INFO("Rough estimates on N50/L50:" << lengths[i - 1] << " / " << i - 1 << " with len " << total_len);
}

PathScore GenomeConsistenceChecker::CountMisassembliesWithStrand(const BidirectionalPath &path, const string strand) const {
    if (strand == "1") {
        return (CountMisassembliesWithStrand(*path.GetConjPath(), "0"));
    }
    PathScore res(0, 0, 0);
    EdgeId prev;
    size_t prev_in_genome = std::numeric_limits<std::size_t>::max();
    size_t prev_in_path = std::numeric_limits<std::size_t>::max();
    MappingRange prev_range;
    for (int i = 0; i < (int) path.Size(); i++) {
        if (genome_spelled_.find(path.At(i)) != genome_spelled_.end()) {
            size_t cur_in_genome =  genome_spelled_[path.At(i)];
            MappingRange cur_range = *gp_.edge_pos.GetEdgePositions(path.At(i), "fxd0").begin();
            if (prev_in_genome != std::numeric_limits<std::size_t>::max()) {
                if (cur_in_genome == prev_in_genome + 1) {
                    int dist_in_genome = (int) cur_range.initial_range.start_pos -  (int) prev_range.initial_range.end_pos;
                    int dist_in_path = (int) path.LengthAt(prev_in_path) - (int) path.LengthAt(i) +  (int) cur_range.mapped_range.start_pos - (int) prev_range.mapped_range.end_pos;
                    DEBUG("Edge " << prev.int_id() << "  position in genome ordering: " << prev_in_genome);
                    DEBUG("Gap in genome / gap in path: " << dist_in_genome << " / " << dist_in_path);
                    if (size_t(abs(dist_in_genome - dist_in_path)) > absolute_max_gap_ && (dist_in_genome * (1 + relative_max_gap_) < dist_in_path || dist_in_path * (1 + relative_max_gap_) < dist_in_genome)) {

                        res.wrong_gap_size ++;
                    }
                } else {
                    if (path.At(i) != circular_edge_ && path.At(prev_in_path) != circular_edge_)
                        res.misassemblies++;
                    else
                        INFO("Skipping fake(circular) misassembly");
                }
            }
            res.mapped_length += cur_range.mapped_range.size();
            prev = path.At(i);
            prev_in_genome = cur_in_genome;
            prev_range = cur_range;
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
    for (auto e: storage_) {
        RefillPos(strand, e);
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
    for(size_t i = 0; i < sz; i++) {
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
    for(size_t i = 0; i < sz; i++) {
        scores[i] = to_process[i].initial_range.size();
        prev[i] = std::numeric_limits<std::size_t>::max();
    }
    for(size_t i = 0; i < sz; i++) {
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
    for(size_t i = 0; i < sz; i++) {
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
    string new_strand = "fxd" + strand;
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
                if (strand == "0") {
                    circular_edge_ = e;
                }
            }
        }

    }
    gp_.edge_pos.AddEdgePosition(e, new_strand, new_mapping);
}



}
