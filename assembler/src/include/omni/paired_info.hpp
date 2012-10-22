//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

//#include "utils.hpp"
//#include "sequence.hpp"
//#include "paired_read.hpp"
#include "io/paired_read.hpp"
#include <cmath>
#include <map>
#include <limits>
#include <xmath.h>
#include "omni_utils.hpp"

//#define MERGE_DATA_RELATIVE_DIFFERENCE 0.3

namespace omnigraph {

/**
 * PairInfo class represents basic data unit for paired information: edges first and second appear
 * in genome at distance d_ and this information has weight weight_.
 */
template<typename EdgeId>
struct PairInfo {
	EdgeId first;
	EdgeId second;
	double d; //distance between starts. Can be negative
	double weight;
	double variance;

	PairInfo(const PairInfo& pair_info) :
			first(pair_info.first), second(pair_info.second), d(pair_info.d), weight(
					pair_info.weight), variance(pair_info.variance) {

	}

	PairInfo(EdgeId first, EdgeId second, double d, double weight,
			double variance) :
			first(first), second(second), d(d), weight(weight), variance(
					variance) {
	}

	/**
	 * Two paired infos are considered equal if they coinside in all parameters except for the weight and the variance.
	 */
	bool operator==(const PairInfo& rhs) const {
		const PairInfo &lhs = *this;
		return lhs.first == rhs.first && lhs.second == rhs.second
				&& math::eq(lhs.d, rhs.d);
	}

	bool operator!=(const PairInfo& rhs) const {
		return !(*this == rhs);
	}

	bool operator<(const PairInfo<EdgeId>& rhs) const {
		const PairInfo<EdgeId>& lhs = *this;
		return lhs.first == rhs.first ?
				lhs.second == rhs.second ?
						math::ls(lhs.d, rhs.d) : lhs.second < rhs.second
				: lhs.first < rhs.first;
	}
};

template<typename Graph>
ostream& operator<<(ostream& os, const PairInfo<Graph>& info) {
	return os << "PairInfo: first=" << info.first << ", second=" << info.second
			<< ", distance=" << info.d << ", weight=" << info.weight
			<< ", variance=" << info.variance;
}

template<typename EdgeId>
const PairInfo<EdgeId> MinPairInfo(EdgeId id) {
	return PairInfo<EdgeId>(id, EdgeId(typename EdgeId::pointer_type(1)),
			-10000000000, 0., 0.);
}

template<typename EdgeId>
const PairInfo<EdgeId> MaxPairInfo(EdgeId id) {
	return PairInfo<EdgeId>(id, EdgeId(typename EdgeId::pointer_type(-1)),
			 10000000000, 0., 0.);
}

template<typename EdgeId>
const PairInfo<EdgeId> MinPairInfo(EdgeId e1, EdgeId e2) {
	PairInfo<EdgeId> info = MinPairInfo(e1);
	info.second = e2;
	return info;
}

template<typename EdgeId>
const PairInfo<EdgeId> MaxPairInfo(EdgeId e1, EdgeId e2) {
	PairInfo<EdgeId> info = MaxPairInfo(e1);
	info.second = e2;
	return info;
}

/**
 * Method returns approximate distance between occurrences of edges in genome rounded to the nearest
 * integer. In case of a tie closest to 0 value is chosen thus one can assume that distance
 * is rounded the same way as opposite one.
 * todo check that written here is true
 */
template<typename EdgeId>
int rounded_d(PairInfo<EdgeId> const& pi) {
	int res = (int) math::round(std::abs(pi.d));
	if (math::ls(pi.d, 0.))
		res = -res;
	return res;
}

template<typename EdgeId>
PairInfo<EdgeId> BackwardInfo(const PairInfo<EdgeId>& pi) {
	return PairInfo<EdgeId>(pi.second, pi.first, -pi.d, pi.weight, pi.variance);
}

template<typename EdgeId>
bool IsSymmetric(PairInfo<EdgeId> const& pi) {
	return pi.first == pi.second && math::eq(pi.d, 0.);
}

//TODO: try storing set<PairInfo>
template <typename EdgeId>
class PairInfoIndexData {
public:
	typedef set<PairInfo<EdgeId> > Data;
	typedef typename Data::iterator data_iterator;
	typedef vector<PairInfo<EdgeId> > PairInfos;

    //  we can not update elements in the std::set, 
    //  although we can delete element, 
    //  and then insert a modified version of it
    //  but here we do not care about safety, and making illegal @const_cast on the std::set element
    void UpdateSingleInfo(const PairInfo<EdgeId>& info, double new_dist, double new_weight, double new_variance) {
        TRACE(info << " is about to be merged with " << new_dist << " " << new_weight << " " << new_variance);
        PairInfo<EdgeId>& info_to_update = const_cast<PairInfo<EdgeId>&>(*data_.find(info));
        using namespace math;
        update_value_if_needed<double>(info_to_update.d, new_dist);
        update_value_if_needed<double>(info_to_update.weight, new_weight);
        update_value_if_needed<double>(info_to_update.variance, new_variance);
    }

    //TODO: rename the method
	void ReplaceFirstEdge(const PairInfo<EdgeId>& info, EdgeId newId) {
		data_.insert(PairInfo<EdgeId>(newId, info.second, info.d, info.weight, info.variance));
	}

	PairInfoIndexData() :
			data_() {
	}

	data_iterator begin() const {
		return data_.begin();
	}

	data_iterator end() const {
		return data_.end();
	}

	size_t size() const {
		return data_.size();
	}

	void AddPairInfo(const PairInfo<EdgeId>& pair_info, bool add_symmetric = 1) {
        data_.insert(pair_info);

		if (!IsSymmetric(pair_info) && add_symmetric)
			data_.insert(BackwardInfo(pair_info));
	}

	void UpdateInfo(const PairInfo<EdgeId>& info, double new_dist, double new_weight, double new_variance, bool add_reversed) {
        // first we update backward info in order to leave @info not modified
		if (add_reversed && !IsSymmetric(info))
			UpdateSingleInfo(BackwardInfo(info), -new_dist, new_weight, new_variance);

		UpdateSingleInfo(info, new_dist, new_weight, new_variance);
	}

	void DeleteEdgeInfo(EdgeId e) {
		set<PairInfo<EdgeId> > paired_edges;

		for (auto lower = LowerBound(e), upper = UpperBound(e); lower != upper;
				++lower) {
			paired_edges.insert(BackwardInfo(*lower));
		}

		for (auto it = paired_edges.begin(); it != paired_edges.end(); ++it) {
			data_.erase(*it);
		}

		data_.erase(LowerBound(e), UpperBound(e));
	}

	void DeletePairInfo(const PairInfo<EdgeId>& info) {
        VERIFY(data_.find(info) != data_.end());
		data_.erase(info);
	}

	void DeleteEdgePairInfo(EdgeId e1, EdgeId e2) {
		data_.erase(LowerBound(e1, e2), UpperBound(e1, e2));
		if (e1 != e2)
			data_.erase(LowerBound(e2, e1), UpperBound(e2, e1));
	}

	PairInfos GetEdgeInfos(EdgeId e) const {
		return PairInfos(LowerBound(e), UpperBound(e));
	}

	PairInfos GetEdgePairInfos(EdgeId e1, EdgeId e2) const {
		return PairInfos(LowerBound(e1, e2), UpperBound(e1, e2));
	}

	void clear() {
		data_.clear();
	}

    data_iterator Find(const PairInfo<EdgeId>& pair_info) const {
        return data_.find(pair_info);   
    }

	data_iterator LowerBound(const PairInfo<EdgeId>& pair_info) const {
		return data_.lower_bound(pair_info);
	}

	data_iterator UpperBound(const PairInfo<EdgeId>& pair_info) const {
		return data_.upper_bound(pair_info);
	}

	data_iterator LowerBound(EdgeId e) const {
		return data_.lower_bound(MinPairInfo(e));
	}

	data_iterator UpperBound(EdgeId e) const {
		return data_.upper_bound(MaxPairInfo(e));
	}

	data_iterator LowerBound(EdgeId e1, EdgeId e2) const {
		return data_.lower_bound(MinPairInfo(e1, e2));
	}

	data_iterator UpperBound(EdgeId e1, EdgeId e2) const {
		return data_.upper_bound(MaxPairInfo(e1, e2));
	}

private:
	Data data_;

    DECL_LOGGER("PairedInfoIndexData");
};

/**
 * PairedInfoIndex stores information about edges connected by paired reads and synchronizes this info with
 * graph.
 */
template<class Graph>
class PairedInfoIndex: public GraphActionHandler<Graph> {

public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef const vector<PairInfo<EdgeId> > PairInfos;
	typedef typename PairInfoIndexData<EdgeId>::data_iterator data_iterator;

public:
	/**
	 * Class EdgePairIterator is used to iterate through paires of edges which have information about distance
	 * between them stored in PairedInfoIndex.
	 */
	class EdgePairIterator {
		typename PairInfoIndexData<EdgeId>::data_iterator position_;
		const PairedInfoIndex<Graph> &index_;
	public:
		EdgePairIterator(
				typename PairInfoIndexData<EdgeId>::data_iterator position,
				const PairedInfoIndex<Graph> &index) :
				position_(position), index_(index) {
		}

		bool operator==(const EdgePairIterator &other) {
			return this->position_ == other.position_;
		}

		bool operator!=(const EdgePairIterator &other) {
			return this->position_ != other.position_;
		}

		PairInfos operator*() const {
			return index_.GetEdgePairInfo(position_->first, position_->second);
		}

		void operator++() {
			position_ = index_.index_data_.UpperBound(position_->first,
					position_->second);
		}

		EdgeId first() const {
		    return position_->first;
		}

		EdgeId second() const {
            return position_->second;
        }
	};

	EdgePairIterator begin() const {
		VERIFY(this->IsAttached());
		return EdgePairIterator(index_data_.begin(), *this);
	}

	EdgePairIterator end() const {
		VERIFY(this->IsAttached());
		return EdgePairIterator(index_data_.end(), *this);
	}

	//begin-end insert size supposed
	PairedInfoIndex(const Graph &g) :
			GraphActionHandler<Graph>(g, "PairedInfoIndex")
    {
	}

	virtual ~PairedInfoIndex() {
		TRACE("~PairedInfoIndex ok");
	}

public:

	void Init() {
		for (auto it = this->g().SmartEdgeBegin(); !it.IsEnd(); ++it) {
			HandleAdd(*it);
		}
	}

	//	void OutputData(ostream &os = cout) {
	//		for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it)
	//			for (auto it1 = graph_.SmartEdgeBegin(); !it1.IsEnd(); ++it1) {
	//				OutputEdgeData(*it, *it1, os);
	//			}
	//	}
	//
	//	void OutputData(string fileName) {
	//		ofstream s;
	//		s.open(fileName.c_str());
	//		OutputData(s);
	//		s.close();
	//	}

	/*
	 * @return quantity of paired info
	 */
	size_t size() const {
		return index_data_.size();
	}

	/**
	 * Method returns all data about given edge
	 */
	PairInfos GetEdgeInfo(EdgeId edge) const {
		VERIFY(this->IsAttached());
		return index_data_.GetEdgeInfos(edge);
	}

	/**
	 * Method returns all data about distance between two given edges
	 */
	PairInfos GetEdgePairInfo(EdgeId first, EdgeId second) const {
		VERIFY(this->IsAttached());
		return index_data_.GetEdgePairInfos(first, second);
	}

	virtual void HandleAdd(EdgeId e) {
		this->AddPairInfo(PairInfo<EdgeId>(e, e, 0., 0., 0.));
	}

	virtual void HandleDelete(EdgeId e) {
		this->RemoveEdgeInfo(e);
	}

	virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
		this->AddPairInfo(PairInfo<EdgeId>(new_edge, new_edge, 0., 0., 0.));
		int shift = 0;
		for (size_t i = 0; i < old_edges.size(); ++i) {
			EdgeId old_edge = old_edges[i];
			TransferInfo(old_edge, new_edge, shift);
			shift -= this->g().length(old_edge);
		}
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		TransferInfo(edge2, new_edge);
		TransferInfo(edge1, new_edge);
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge1,
			EdgeId new_edge2) {
		double prop = (double) this->g().length(new_edge1)
				/ this->g().length(old_edge);
		//		size_t shift = graph_.length(new_edge1);
		TransferInfo(old_edge, new_edge1, 0, prop);
		//		PassEdge(graph_.length(new_edge1), shift);
		TransferInfo(old_edge, new_edge2, this->g().length(new_edge1),
				1 - prop);
	}

	//	bool Check() {
	//		for (auto it = index_data_.begin(); it != index_data_.end(); ++it) {
	//			PairInfo<EdgeId> inf = *it;
	//			if (inf.first != inf.second && inf.d > 0) {
	////				cout << inf.first << endl;
	////				cout << graph_.length(inf.first) << endl;
	//				if (inf.d <= graph_.length(inf.first) - 1) {
	//					cout << "gopa " << inf.first << " " << inf.second << endl;
	//					cout << graph_.length(inf.first) << " " << graph_.length(inf.second) << endl;
	//					cout << graph_.EdgeStart(inf.first) << " " << graph_.EdgeEnd(inf.first);
	//					cout << graph_.EdgeStart(inf.second) << " " << graph_.EdgeEnd(inf.second);
	//					VERIFY(false);
	//					return false;
	//				}
	//			}
	//		}
	//		return true;
	//	}

	/**
	 * Method allows to add pair info to index directly instead of filling it from stream.
     * Notice, that you can merge two pair infos only if their variances are equal to zero!
	 */
	void AddPairInfo(const PairInfo<EdgeId>& pair_info, bool add_reversed = 1) {
		VERIFY(this->IsAttached());
		TRACE("Adding pair info to pi index: " << pair_info.first << pair_info.second << " " << index_data_.size());

		data_iterator cluster_it = index_data_.Find(pair_info);
        if (cluster_it != index_data_.end()) {
            TRACE("Such pair info exists, merging now");
            const PairInfo<EdgeId>& existing_info = *cluster_it;
            VERIFY(existing_info == pair_info);
			MergeData(existing_info, pair_info, add_reversed);
		} else {
            TRACE("Such pair info does not exist");
			index_data_.AddPairInfo(pair_info, add_reversed);
		}
	}

    void AddAll(const PairedInfoIndex<Graph>& paired_index) {
		VERIFY(this->IsAttached());
        for (auto iter = paired_index.begin(); iter != paired_index.end(); ++iter) {
            const vector<PairInfo<EdgeId> >& infos = *iter;
            for (auto pi_iter = infos.begin(); pi_iter != infos.end(); ++pi_iter) {
                this->AddPairInfo(*pi_iter, false);
            }
        }
    }

	void RemoveEdgeInfo(EdgeId edge) {
		VERIFY(this->IsAttached());
		index_data_.DeleteEdgeInfo(edge);
	}

	void RemovePairInfo(const PairInfo<EdgeId>& pair_info) {
		VERIFY(this->IsAttached());
		index_data_.DeletePairInfo(pair_info);
	}

	void Clear() {
	    index_data_.clear();
	}

private:
	PairInfoIndexData<EdgeId> index_data_;

	void MergeData(const PairInfo<EdgeId>& info_to_update, const PairInfo<EdgeId>& info_to_add,
			bool add_reversed) {
        double weight_to_add = info_to_add.weight;

        // counting new bounds in the case, when we are merging pair infos with non-zero variance
        double left_bound = min(info_to_update.d - info_to_update.variance, 
                info_to_add.d - info_to_add.variance);
        double right_bound = max(info_to_update.d + info_to_update.variance, 
                info_to_add.d + info_to_add.variance);
        double new_dist = (left_bound + right_bound) * 0.5;
        double new_weight = info_to_update.weight + weight_to_add;
        double new_variance = (right_bound - left_bound) * 0.5;

        index_data_.UpdateInfo(info_to_update, new_dist, new_weight, new_variance, add_reversed);
	}

	//	void OutputEdgeData(EdgeId edge1, EdgeId edge2, ostream &os = cout) {
	//		PairInfos vec = GetEdgePairInfo(edge1, edge2);
	//		if (vec.size() != 0) {
	//			os << edge1 << " " << this->g().length(edge1) << " " << edge2 << " "
	//					<< this->g().length(edge2) << endl;
	//			if (this->g().EdgeEnd(edge1) == this->g().EdgeStart(edge2))
	//				os << "+" << endl;
	//			if (this->g().EdgeEnd(edge2) == this->g().EdgeStart(edge1))
	//				os << "-" << endl;
	//			int min = INT_MIN;
	//			for (size_t i = 0; i < vec.size(); i++) {
	//				int next = -1;
	//				for (size_t j = 0; j < vec.size(); j++) {
	//					if (vec[j].d > min
	//							&& (next == -1 || vec[next].d > vec[j].d)) {
	//						next = j;
	//					}
	//				}
	//				os << vec[next].d << " " << vec[next].weight << endl;
	//				if (next == -1) {
	//					VERIFY(false);
	//				}
	//				if (vec[next].d > 100000) {
	//					VERIFY(false);
	//				}
	//				min = vec[next].d;
	//			}
	//		}
	//	}

	void TransferInfo(EdgeId old_edge, EdgeId new_edge, int shift = 0,
			double weight_scale = 1.0) {
		PairInfos pair_infos = GetEdgeInfo(old_edge);
		for (size_t j = 0; j < pair_infos.size(); ++j) {
			PairInfo<EdgeId> old_pair_info = pair_infos[j];
			if (old_edge != old_pair_info.second) {
				AddPairInfo(
						PairInfo<EdgeId>(new_edge, old_pair_info.second,
								old_pair_info.d - shift,
								weight_scale * old_pair_info.weight,
								old_pair_info.variance));
			} else if (!math::eq(old_pair_info.d, 0.)) {
				AddPairInfo(
						PairInfo<EdgeId>(new_edge, new_edge, old_pair_info.d,
								weight_scale * 0.5 * old_pair_info.weight,
								old_pair_info.variance));
			} else {
				AddPairInfo(
						PairInfo<EdgeId>(new_edge, new_edge, old_pair_info.d,
								weight_scale * old_pair_info.weight,
								old_pair_info.variance));
			}
		}
	}

	size_t CorrectLength(const Path<EdgeId>& path, size_t index) const {
		size_t result = this->g().length(path[index]);
		if (index == 0) {
			result -= path.start_pos();
		}
		if (index == path.size() - 1) {
			result -= this->g().length(path[index]) - path.end_pos();
		}
		return result;
	}

    DECL_LOGGER("PairedInfoIndex");
};

template<class Graph, class SequenceMapper, class PairedStream>
class PairedIndexFiller {
private:
	typedef typename Graph::EdgeId EdgeId;

	const Graph &graph_;

	const SequenceMapper& mapper_;

	std::vector< PairedStream* > streams_;

	size_t CorrectLength(Path<EdgeId> path, size_t idx) {
		size_t answer = graph_.length(path[idx]);
		if (idx == 0)
			answer -= path.start_pos();
		if (idx == path.size() - 1)
			answer -= graph_.length(path[idx]) - path.end_pos();
		return answer;
	}

	template<class PairedRead>
	void ProcessPairedRead(omnigraph::PairedInfoIndex<Graph> &paired_index,
			const PairedRead& p_r) {
		Sequence read1 = p_r.first().sequence();
		Sequence read2 = p_r.second().sequence();
		Path<EdgeId> path1 = mapper_.MapSequence(read1);
		Path<EdgeId> path2 = mapper_.MapSequence(read2);
		size_t distance = p_r.distance();
		int current_distance1 = distance + path1.start_pos()
				- path2.start_pos();
		for (size_t i = 0; i < path1.size(); ++i) {
			int current_distance2 = current_distance1;
			for (size_t j = 0; j < path2.size(); ++j) {
				double weight = CorrectLength(path1, i)
						* CorrectLength(path2, j);
				PairInfo<EdgeId> new_info(path1[i], path2[j], current_distance2,
						weight, 0.);
				paired_index.AddPairInfo(new_info);
				current_distance2 += graph_.length(path2[j]);
			}
			current_distance1 -= graph_.length(path1[i]);
		}
	}

    /**
     * Method reads paired data from stream, maps it to genome and stores it in this PairInfoIndex.
     */
    void FillUsualIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
        for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            paired_index.AddPairInfo(PairInfo<EdgeId>(*it, *it, 0., 0., 0.));
        }

        INFO("Processing paired reads (takes a while)");

        PairedStream& stream = *(streams_.front());
        stream.reset();
        size_t n = 0;
        while (!stream.eof()) {
            typename PairedStream::read_type p_r;
            stream >> p_r;
            ProcessPairedRead(paired_index, p_r);
            VERBOSE_POWER(++n, " paired reads processed");
        }
    }


    void FillParallelIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
        for (auto it = graph_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            paired_index.AddPairInfo(PairInfo<EdgeId>(*it, *it, 0, 0.0, 0.));
        }

        INFO("Processing paired reads (takes a while)");

        size_t nthreads = streams_.size();
        std::vector< omnigraph::PairedInfoIndex<Graph>* > buffer_pi(nthreads);
        buffer_pi[0] = &paired_index;

        for (size_t i = 1; i < nthreads; ++i) {
            buffer_pi[i] = new omnigraph::PairedInfoIndex<Graph>(graph_, paired_index.GetMaxDifference());
        }

        size_t counter = 0;
        #pragma omp parallel num_threads(nthreads)
        {
            #pragma omp for reduction(+ : counter)
            for (size_t i = 0; i < nthreads; ++i) {

                typename PairedStream::read_type r;
                PairedStream& stream = *streams_[i];
                stream.reset();

                while (!stream.eof()) {
                    stream >> r;
                    ++counter;

                    ProcessPairedRead(*buffer_pi[i], r);
                }
            }
        }

        INFO("Used " << counter << " paired reads");

        INFO("Merging paired indices");
        for (size_t i = 1; i < nthreads; ++i) {
            buffer_pi[0]->AddAll(*(buffer_pi[i]));
            delete buffer_pi[i];
        }

    }

public:

	PairedIndexFiller(const Graph &graph, const SequenceMapper& mapper,
	        PairedStream& stream) :
			graph_(graph), mapper_(mapper), streams_() {

	    streams_.push_back(&stream);
	}

    PairedIndexFiller(const Graph &graph, const SequenceMapper& mapper,
            std::vector<PairedStream*>& streams) :
            graph_(graph), mapper_(mapper), streams_(streams.begin(), streams.end()) {
    }

    void FillIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
        if (streams_.size() == 1) {
            FillUsualIndex(paired_index);
        } else {
            FillParallelIndex(paired_index);
        }
    }


private:
	DECL_LOGGER("PairedIndexFiller");

};

//New metric weight normalizer
template<class Graph>
class PairedInfoWeightNormalizer {
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	const size_t insert_size_;
	//todo use this param!
	const double is_var_;
	const size_t read_length_;
	const size_t k_;
	const double avg_coverage_;
public:

	//Delta better to be around 5-10% of insert size
	PairedInfoWeightNormalizer(const Graph& g, size_t insert_size,
			double is_var, size_t read_length, size_t k, double avg_coverage) :
			g_(g), insert_size_(insert_size), is_var_(is_var), read_length_(
					read_length), k_(k), avg_coverage_(avg_coverage) {
	}

    const PairInfo<EdgeId> NormalizeWeightWithCoverage(const PairInfo<EdgeId> & pair_info) {
        PairInfo<EdgeId> new_info = pair_info;
        new_info.weight *= g_.length(pair_info.first) * g_.length(pair_info.second) * 1. / (g_.coverage(pair_info.first) * g_.coverage(pair_info.second)); 
        return new_info;
    }

	const PairInfo<EdgeId> NormalizeWeight(const PairInfo<EdgeId>& pair_info) {
		double w = 0.;
		if (math::eq(pair_info.d, 0.) && pair_info.first == pair_info.second) {
			w = 0. + g_.length(pair_info.first) - insert_size_
					+ 2. * read_length_ + 1. - k_;
		} else {
			EdgeId e1 =
					(math::ge(pair_info.d, 0.)) ?
							pair_info.first : pair_info.second;
			EdgeId e2 =
					(math::ge(pair_info.d, 0.)) ?
							pair_info.second : pair_info.first;
			int gap_len = std::abs(rounded_d(pair_info)) - g_.length(e1);
			int right = std::min(insert_size_,
					gap_len + g_.length(e2) + read_length_);
			int left = std::max(gap_len,
					int(insert_size_) - int(read_length_) - int(g_.length(e1)));
			w = 0. + right - left + 1 - k_;
		}

		double result_weight = pair_info.weight;
		if (math::gr(w, /*-10.*/0.)) {
			result_weight /= w; //(w + 10);
		} else result_weight = 0;
		double cov_norm_coeff = avg_coverage_ / (2. * (read_length_ - k_));
		result_weight /= cov_norm_coeff;

		PairInfo<EdgeId> result(pair_info);
		result.weight = result_weight;
		return result;
	}
};

template<class Graph>
class JumpingNormilizerFunction {
private:
	typedef typename Graph::EdgeId EdgeId;
	const Graph& graph_;
	size_t read_length_;
	size_t max_norm_;

public:
	JumpingNormilizerFunction(const Graph& graph, size_t read_length,
			size_t max_norm) :
			graph_(graph), read_length_(read_length), max_norm_(max_norm) {
	}

	size_t norm(EdgeId first, EdgeId second) const {
		return std::min(std::min(graph_.length(first), graph_.length(second)),
				max_norm_) + read_length_ - graph_.k();
	}

	const PairInfo<EdgeId> operator()(const PairInfo<EdgeId>& pair_info) const {
		return PairInfo<EdgeId>(pair_info.first, pair_info.second, pair_info.d,
				pair_info.weight / norm(pair_info.first, pair_info.second),
				pair_info.variance);
	}
};

template<class Graph>
const PairInfo<typename Graph::EdgeId> TrivialWeightNormalization(
		const PairInfo<typename Graph::EdgeId>& pair_info) {
	return pair_info;
}

template<class Graph>
class PairedInfoNormalizer {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef boost::function<const PairInfo<EdgeId>(const PairInfo<EdgeId>&)> WeightNormalizer;
private:

	WeightNormalizer normalizing_function_;
public:

	PairedInfoNormalizer(WeightNormalizer normalizing_function) :
			normalizing_function_(normalizing_function) {

	}

// temporary due to path_extend absolute thresholds
	void FillNormalizedIndex(const PairedInfoIndex<Graph>& paired_index, PairedInfoIndex<Graph>& normalized_index, double coeff = 1.) const {
		for (auto it = paired_index.begin(); it != paired_index.end(); ++it) {
			const vector<PairInfo<EdgeId> >& infos = *it;
			for (auto it2 = infos.begin(); it2 != infos.end(); ++it2) {
                PairInfo<EdgeId> tmp = *it2;
                tmp.weight *= coeff;
                DEBUG("Normalized pair info " << (tmp) << " " << normalizing_function_(tmp));
				normalized_index.AddPairInfo(normalizing_function_(tmp), false);
			}
		}
	}
};

//template<class Graph>
//class PairedInfoSymmetryHack {
//public:
//	typedef typename Graph::EdgeId EdgeId;
//private:
//	const Graph& graph_;
//	const PairedInfoIndex<Graph>& paired_index_;
//public:
//
//	PairedInfoSymmetryHack(const Graph& graph,
//			const PairedInfoIndex<Graph>& paired_index) :
//			graph_(graph), paired_index_(paired_index) {
//	}
//
//	void FillSymmetricIndex(PairedInfoIndex<Graph>& index) {
//		for (auto it = paired_index_.begin(); it != paired_index_.end(); ++it) {
//			vector<PairInfo<EdgeId>> infos = *it;
//			for (auto it2 = infos.begin(); it2 != infos.end(); ++it2) {
//				index.AddPairInfo(
//						PairInfo<EdgeId>(it2->first, it2->second, it2->d,
//								it2->weight * 0.5, it2->variance), 0);
//				index.AddPairInfo(
//						PairInfo<EdgeId>(
//								graph_.conjugate(it2->second),
//								graph_.conjugate(it2->first),
//								it2->d + graph_.length(it2->second)
//										- graph_.length(it2->first),
//								it2->weight * 0.5, it2->variance), 0);
//				auto info = index.GetEdgePairInfo(infos[0].first,
//						infos[0].second);
//				auto symmetric_info = index.GetEdgePairInfo(
//						graph_.conjugate(info[0].second),
//						graph_.conjugate(info[0].first));
//			}
//		}
//	}
//};

}
