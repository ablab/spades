/*
 * pac_index.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: lab42
 */
#pragma once

#include "indices/perfect_hash_map.hpp"
#include "graph_pack.hpp"
#include <algorithm>
using std::map;
using std::set;
namespace pacbio {
template<class T>
struct pair_iterator_less {
	bool operator ()(pair<size_t, T> const& a, pair<size_t, T> const& b) const {
		return (a.first < b.first);
	}
};

struct MappingInstance {
	int edge_position;
	int read_position;
	//Now quality is the same with multiplicity, so best quality is 1,
	int quality;
	MappingInstance(int edge_position, int read_position, int quality) :
			edge_position(edge_position), read_position(read_position), quality(quality) {
	}

	inline bool IsUnique() const {
		return (quality == 1);
	}

	string str() {
		stringstream s;
		s << "E: " << edge_position << " R: " << read_position << " Q: " << quality;
		return s.str();
	}

//Less by EDGE position
	bool operator <(MappingInstance const& b) const {
		if (edge_position < b.edge_position || (edge_position == b.edge_position && read_position < b.read_position))
			return true;
		else
			return false;
	}
};

//Less by READ position
struct ReadPositionComparator {
	bool operator ()(MappingInstance const& a, MappingInstance const& b) const {
		return (a.read_position < b.read_position || (a.read_position == b.read_position && a.edge_position < b.edge_position));
	}
};

template<class Graph>
struct KmerCluster {
	typedef typename Graph::EdgeId EdgeId;
	int last_trustable_index;
	int first_trustable_index;
	EdgeId edgeId;
	vector<MappingInstance> sorted_positions;
	int size;

	KmerCluster(EdgeId e, const vector<MappingInstance>& v) {
		last_trustable_index = 0;
		first_trustable_index = 0;
		edgeId = e;
		size = (int) v.size();
		sorted_positions = v;
		FillTrustableIndeces();
	}

	bool operator <(const KmerCluster & b) const {
		return (edgeId < b.edgeId || (edgeId == b.edgeId && sorted_positions < b.sorted_positions));
	}

	bool CanFollow(const KmerCluster &b) const {
		return (b.sorted_positions[b.last_trustable_index].read_position < sorted_positions[first_trustable_index].read_position);
	}

	void FillTrustableIndeces() {
		//ignore non-unique kmers for distance determination
		int first_unique_ind = 0;
		while (first_unique_ind != size - 1 && !(sorted_positions[first_unique_ind].IsUnique())) {
			first_unique_ind += 1;
		}
		int last_unique_ind = size - 1;
		while (last_unique_ind != 0 && !(sorted_positions[last_unique_ind].IsUnique())) {
			last_unique_ind -= 1;
		}
		last_trustable_index = last_unique_ind;
		first_trustable_index = first_unique_ind;
	}

    string str(const Graph &g) const{
        stringstream s;
        s << "Edge: " << g.int_id(edgeId) << " " << first_trustable_index<< " - "  <<last_trustable_index << size;
        return s.str();
    }
};

template<class Graph>
struct GapDescription {
	typedef typename Graph::EdgeId EdgeId;
	typename Graph::EdgeId start, end;
	Sequence gap_seq;
	int edge_gap_start_position, edge_gap_end_position;


	GapDescription(EdgeId start_e, EdgeId end_e, const Sequence &gap, int gap_start, int gap_end) :
			start(start_e), end(end_e), gap_seq(gap.str()), edge_gap_start_position(gap_start), edge_gap_end_position(gap_end) {
	}

	GapDescription(const KmerCluster<Graph> &a, const KmerCluster<Graph> & b, Sequence read, int pacbio_k) {
		edge_gap_start_position = a.sorted_positions[a.last_trustable_index].edge_position;
		edge_gap_end_position = b.sorted_positions[b.first_trustable_index].edge_position + pacbio_k - 1;
		start = a.edgeId;
		end = b.edgeId;
		DEBUG(read.str());
		gap_seq = read.Subseq(a.sorted_positions[a.last_trustable_index].read_position, b.sorted_positions[b.first_trustable_index].read_position + pacbio_k - 1);
		DEBUG(gap_seq.str());
		DEBUG("gap added");
	}

	GapDescription<Graph> conjugate(Graph &g_, int shift) const {
        GapDescription<Graph> res(
                g_.conjugate(end), g_.conjugate(start), (!gap_seq),
                (int) g_.length(end) + shift - edge_gap_end_position,
                (int) g_.length(start) + shift - edge_gap_start_position);
		 DEBUG("conjugate created" << res.str(g_));
		 return res;
	}

	string str(Graph &g_) const {
		stringstream s;
		s << g_.int_id(start) << " " << edge_gap_start_position <<endl << g_.int_id(end) << " " << edge_gap_end_position << endl << gap_seq.str()<< endl;
		return s.str();
	}

	bool operator <(const GapDescription & b) const {
		return (start < b.start || (start == b.start &&  end < b.end) ||
				(start == b.start &&  end == b.end && edge_gap_start_position < b.edge_gap_start_position));
	}

private:
	DECL_LOGGER("PacIndex")
	;
};

template<class Graph>
struct OneReadMapping {
	typedef typename Graph::EdgeId EdgeId;
	vector<vector<EdgeId> > main_storage;
	vector<GapDescription<Graph> > gaps;
	vector<size_t> real_length;
	OneReadMapping(vector<vector<EdgeId> > &paths_description, vector<GapDescription<Graph> > &gaps_description, vector<size_t> real_length) :
			main_storage(paths_description), gaps(gaps_description), real_length(real_length) {
	}

};


struct StatsCounter{

    map<size_t,size_t> path_len_in_edges;
    vector<size_t> subreads_length;
    size_t total_len ;
    size_t reads_with_conjugate;
    size_t subreads_count;
    StatsCounter() {
        total_len = 0;
        reads_with_conjugate = 0;
    }

    void AddStorage(StatsCounter &other) {
        total_len += other.total_len;
        reads_with_conjugate += other.reads_with_conjugate;
        for (auto iter = other.subreads_length.begin(); iter != other.subreads_length.end(); ++iter) {
            subreads_length.push_back(*iter);
        }
        for (auto iter = other.path_len_in_edges.begin(); iter != other.path_len_in_edges.end(); ++iter){
            auto j_iter = iter;
            if (( j_iter = path_len_in_edges.find(iter->first)) == other.path_len_in_edges.end()){
                path_len_in_edges.insert(make_pair(iter->first, iter->second));
            } else {
                path_len_in_edges[j_iter->first] += iter->second;
            }
        }
    }

};

inline int StringDistance(string &a, string &b) {
    int a_len = (int) a.length();
    int b_len = (int) b.length();
    int d = min(a_len / 3, b_len / 3);
    d = max(d, 10);
    DEBUG(a_len << " " << b_len << " " << d);
    vector<vector<int> > table(a_len);
    //int d =
    for (int i = 0; i < a_len; i++) {
        table[i].resize(b_len);
        int low = max(max(0, i - d - 1), i + b_len - a_len - d - 1);
        int high = min(min(b_len, i + d + 1), i + a_len - b_len + d + 1);
        TRACE(low << " " <<high);
        for (int j = low; j < high; j++)
            table[i][j] = 1000000;
    }
    table[a_len - 1][b_len - 1] = 1000000;
    table[0][0] = 0;
//free deletions on begin
//      for(int j = 0; j < b_len; j++)
//          table[0][j] = 0;

    for (int i = 0; i < a_len; i++) {
        int low = max(max(0, i - d), i + b_len - a_len - d);
        int high = min(min(b_len, i + d), i + a_len - b_len + d);

        TRACE(low << " " <<high);
        for (int j = low; j < high; j++) {

            if (i > 0)
                table[i][j] = min(table[i][j], table[i - 1][j] + 1);
            if (j > 0)
                table[i][j] = min(table[i][j], table[i][j - 1] + 1);
            if (i > 0 && j > 0) {
                int add = 1;
                if (a[i] == b[j])
                    add = 0;
                table[i][j] = min(table[i][j], table[i - 1][j - 1] + add);
            }
        }
    }
    //return table[a_len - 1][b_len - 1];
//free deletions on end
    int res = table[a_len - 1][b_len - 1];
    DEBUG(res);
//      for(int j = 0; j < b_len; j++){
//          res = min(table[a_len - 1][j], res);
//      }
    return res;
}


}
