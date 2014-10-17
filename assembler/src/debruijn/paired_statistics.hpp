#pragma once

template<class Graph>
class PairInfoChecker {
private:
	typedef typename Graph::EdgeId EdgeId;
	const EdgesPositionHandler<Graph> &positions_;
	size_t first_bound_;
	const size_t second_bound_;
	vector<double> perfect_matches_;
	vector<double> good_matches_;
	vector<double> mismatches_;
	vector<double> imperfect_matches_;

public:
	PairInfoChecker(const EdgesPositionHandler<Graph> &positions,
			size_t first_bound, size_t second_bound) :
			positions_(positions), first_bound_(first_bound), second_bound_(
					second_bound) {
	}

	void Check(const de::PairedInfoIndex<Graph> &paired_index) {
		for (auto it = paired_index.begin(); it != paired_index.end(); ++it) {
			auto vec = *it;
			for (auto vec_it = vec.begin(); vec_it != vec.end(); ++vec_it) {
				size_t code = CheckSingleInfo(*vec_it);
				if (code == 0) {
					perfect_matches_.push_back(vec_it->weight);
				} else if (code == 1) {
					good_matches_.push_back(vec_it->weight);
				} else if (code == 2) {
					mismatches_.push_back(vec_it->weight);
				} else if (code == 3) {
					imperfect_matches_.push_back(vec_it->weight);
				}
			}
		}
	}

	size_t CheckSingleInfo(de::PairInfo<EdgeId> info) {
		const vector<EdgePosition> &pos1 = positions_.GetEdgePositions(
				info.first);
		const vector<EdgePosition> &pos2 = positions_.GetEdgePositions(
				info.second);
		bool good_match_found = false;
		for (size_t i = 0; i < pos1.size(); i++)
			for (size_t j = 0; j < pos2.size(); j++) {
				if (abs(pos1[i].mr.initial_range.start_pos + info.d - pos2[j].mr.initial_range.start_pos)
						<= first_bound_ + info.variance) {
					if (info.variance == 0) {
						return 0;
					} else {
						return 3;
					}
				} else if (abs(pos1[i].mr.initial_range.start_pos + info.d - pos2[j].mr.initial_range.start_pos)
						<= second_bound_) {
					good_match_found = true;
				}
			}
		if (good_match_found) {
			return 1;
		} else {
			return 2;
		}
	}

	void WriteResultsToFile(vector<double> results, const string &file_name) {
		sort(results.begin(), results.end());
		ofstream os;
		os.open(file_name.c_str());
		for (size_t i = 0; i < results.size(); i++) {
			os << results[i] << endl;
		}
		os.close();
	}

	void WriteResults(const string &folder_name) {
        path::make_dir(folder_name);
		WriteResultsToFile(perfect_matches_,
				folder_name + "/perfect_matches.txt");
		WriteResultsToFile(good_matches_, folder_name + "/good_matches.txt");
		WriteResultsToFile(mismatches_, folder_name + "/mismatches.txt");
		WriteResultsToFile(imperfect_matches_,
				folder_name + "/imperfect_matches.txt");
	}
};

template<class Graph>
class TrivialEdgePairChecker {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    const size_t bound_;
public:
    TrivialEdgePairChecker(const Graph &graph, size_t bound = (size_t) - 1) :
            graph_(graph), bound_(bound) {
    }

    /*
     * Very bad code. Shame on me.
     */
    bool GoForward(EdgeId &edge) {
        if (!graph_.CheckUniqueOutgoingEdge(graph_.EdgeEnd(edge))) {
            return false;
        }
        edge = graph_.GetUniqueOutgoingEdge(graph_.EdgeEnd(edge));
        return true;
    }

    bool GoBackward(EdgeId &edge) {
        if (!graph_.CheckUniqueIncomingEdge(graph_.EdgeStart(edge))) {
            return false;
        }
        edge = graph_.GetUniqueIncomingEdge(graph_.EdgeStart(edge));
        return true;
    }

    bool CheckForward(EdgeId edge1, EdgeId edge2) {
        set<EdgeId> was;
        size_t length = 0;
        do {
            if (edge1 == edge2)
                return true;
            if (was.count(edge1) != 0)
                return false;
            was.insert(edge1);
            length += graph_.length(edge1);
        } while (length <= bound_ && GoForward(edge1));
        return false;
    }

    bool CheckBackward(EdgeId edge1, EdgeId edge2) {
        set<EdgeId> was;
        size_t length = 0;
        do {
            if (edge1 == edge2)
                return true;
            if (was.count(edge1) != 0)
                return false;
            was.insert(edge1);
            length += graph_.length(edge1);
        } while (length <= bound_ && GoBackward(edge1));
        return false;
    }

    bool Check(EdgeId edge1, EdgeId edge2) {
        return CheckForward(edge1, edge2) || CheckBackward(edge2, edge1)
        /*|| CheckForward(edge2, edge1) || CheckBackward(edge1, edge2)*/;
    }
};

template<class Graph>
class EdgePairStat: public AbstractStatCounter {

private:
  typedef typename Graph::EdgeId EdgeId;
  typedef pair<EdgeId, EdgeId> EdgePair;
  const Graph& graph_;
  const PairedInfoIndexT<Graph>& pair_info_;
  const string& output_folder_;

public:
  EdgePairStat(const Graph &graph, const PairedInfoIndexT<Graph> &pair_info,
      const string &output_folder) :
      graph_(graph), pair_info_(pair_info), output_folder_(output_folder) {
  }

  virtual ~EdgePairStat() {
  }

private:
  vector<double> GetWeights(map<EdgePair, double>& edge_pairs) {
    vector<double> weights;
    for (auto it = edge_pairs.begin(); it != edge_pairs.end(); ++it) {
      weights.push_back(it->second);
    }
    sort(weights.begin(), weights.end());
    return weights;
  }

  void GetPairInfo(map<EdgePair, double> &edge_pairs, PairedInfoIndexT<Graph>& index) {
    for (auto it = index.begin(); it != index.end(); ++it) {
      de::Histogram v = *it;
      size_t w = 0;
      for (auto I = v.begin(); I != v.end(); ++I)
        w += (size_t) I->weight;

      edge_pairs.insert(make_pair(make_pair(it.first(), it.second()), w));
    }
  }

  void RemoveTrivial(map<pair<EdgeId, EdgeId> , double> &edge_pairs) {
    TrivialEdgePairChecker<Graph> checker(graph_);
    for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end();
        ) {
      if (checker.Check(iterator->first.first, iterator->first.second)) {
        edge_pairs.erase(iterator++);
      } else {
        ++iterator;
      }
    }
  }

  //  void RemoveUntrustful(map<pair<EdgeId, EdgeId> , double> &edge_pairs, double bound) {
  //    vector<double> weights;
  //    for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end(); ++iterator) {
  //      weights.push_back(iterator->second);
  //    }
  //    sort(weights.begin(), weights.end());
  //
  //    for (auto iterator = edge_pairs.begin(); iterator != edge_pairs.end();) {
  //      if(iterator->second < bound) {
  //        edge_pairs.erase(iterator++);
  //      } else {
  //        ++iterator;
  //      }
  //    }
  //  }

public:
  vector<pair<int, double>> ComulativeHistogram(vector<double> weights) {
    vector<pair<int, double>> result;
    int cur = weights.size() - 1;
    size_t max = 1000;
    vector<double> res(max);
    for (int i = max - 1; i >= 0; i--) {
      while (cur >= 0 && weights[cur] >= i + 1) {
        cur--;
      }
      res[i] = weights.size() - 1 - cur;
    }
    for (size_t i = 0; i < weights.size(); i++) {
      result.push_back(make_pair(i + 1, res[i]));
    }
    return result;
  }

  //  void OutputWeights(vector<double> weights, string file_name) {
  //    ofstream os(file_name);
  //    size_t cur = weights.size() - 1;
  //    size_t max = 1000;
  //    vector<double> res(max);
  //    for(int i = max - 1; i >= 0; i--) {
  //      while(cur >= 0 && weights[cur] >= i + 1) {
  //        cur--;
  //      }
  //      res[i] = weights.size() - 1 - cur;
  //    }
  //    for(size_t i = 0; i < weights.size(); i++) {
  //      os << i + 1 << " " << res[i] << endl;
  //    }
  //    os.close();
  //  }

  bool ContainsPositiveDistance(EdgeId e1, const Histogram& infos) const {
    int first_len = int(graph_.length(e1));
    for (auto it = infos.begin(); it != infos.end(); ++it) {
      if (rounded_d(*it) > first_len)
        return true;
    }
    return false;
  }

  virtual void Count() {
    typedef pair<EdgeId, EdgeId> EdgePair;
    PairedInfoIndexT<Graph> new_index = pair_info_;
    PairInfoWeightFilter<Graph>(graph_, 40).Filter(new_index);
    map<EdgePair, double> edge_pairs;
    TrivialEdgePairChecker<Graph> checker(graph_);
    size_t nontrivial = 0;
    size_t pair_number = 0;
    for (auto iterator = new_index.begin(); iterator != new_index.end(); ++iterator) {
      Histogram info = *iterator;
      if (ContainsPositiveDistance(iterator.first(), info)) {
        ++pair_number;
        if (checker.Check(iterator.first(), iterator.second())) {
          ++nontrivial;
        }
      }
    }
    GetPairInfo(edge_pairs, new_index);
    INFO("Number of edge pairs connected with paired info: " << pair_number);
    RemoveTrivial(edge_pairs);
    INFO("Number of nontrivial edge pairs connected with paired info: " << nontrivial);
  }
};

template<class Graph>
class UniquePathStat: public AbstractStatCounter {

  typedef typename Graph::EdgeId EdgeId;
  const Graph& g_;
  const PairedInfoIndexT<Graph>& filtered_index_;
  size_t insert_size_;
  size_t max_read_length_;
  size_t gap_;
  double variance_delta_;

  size_t considered_edge_pair_cnt_;
  size_t unique_distance_cnt_;
  size_t non_unique_distance_cnt_;

  bool ContainsPositiveDistance(EdgeId e1, const de::Histogram& infos) const {
    int first_len = int(g_.length(e1));
    for (auto it = infos.begin(); it != infos.end(); ++it) {
      if (rounded_d(*it) > first_len)
        return true;
    }
    return false;
  }

public:

  UniquePathStat(const Graph& g, const PairedInfoIndexT<Graph>& filtered_index,
      size_t insert_size, size_t max_read_length, double variance_delta) :
      g_(g), filtered_index_(filtered_index), insert_size_(insert_size), max_read_length_(
          max_read_length), gap_(insert_size_ - 2 * max_read_length_), variance_delta_(
          variance_delta), considered_edge_pair_cnt_(0), unique_distance_cnt_(
          0), non_unique_distance_cnt_(0) {

  }

  virtual void Count() {
    //    PairedInfoIndexT<Graph> filtered_index(g_);
    //    PairInfoFilter < Graph > (g_, 40).Filter(pair_info_, filtered_index);

    for (auto it = filtered_index_.begin(); it != filtered_index_.end(); ++it)
    {
      if (ContainsPositiveDistance(it.first(), *it)) {
        considered_edge_pair_cnt_++;
        EdgeId e1 = it.first();
        EdgeId e2 = it.second();

        //        cout << "Finding paths between edges " << e1 << " and " << e2 << endl;
        NonEmptyPathCounter<Graph> counter(g_);
        //        VertexLablerCallback<Graph> graph_labeler(g_);
        //        CompositeCallback<Graph> composite_callback;
        //        composite_callback.AddProcessor(counter);
        //        composite_callback.AddProcessor(graph_labeler);
        PathProcessor<Graph> path_processor(
            g_,
            omnigraph::PairInfoPathLengthLowerBound(g_.k(),
                g_.length(e1), g_.length(e2), (int) gap_,
                variance_delta_),
            omnigraph::PairInfoPathLengthUpperBound(g_.k(),
                insert_size_, variance_delta_),
            g_.EdgeEnd(e1),
            g_.EdgeStart(e2), 
            counter);
        path_processor.Process();
        if (counter.count() == 1) {
          unique_distance_cnt_++;
        }
        if (counter.count() > 1) {
          non_unique_distance_cnt_++;

        }
      }
    }
    INFO("Considered " << considered_edge_pair_cnt_ << " edge pairs")INFO(
        unique_distance_cnt_ << " edge pairs connected with unique path of appropriate length")
    INFO(
        non_unique_distance_cnt_ << " edge pairs connected with non-unique path of appropriate length")
  }

  size_t considered_edge_pair_count() {
    return considered_edge_pair_cnt_;
  }

  size_t unique_distance_count() {
    return unique_distance_cnt_;
  }

  size_t non_unique_distance_count() {
    return non_unique_distance_cnt_;
  }
private:
  DECL_LOGGER("UniquePathStat")
};

template<class Graph>
class MatePairTransformStat: public AbstractStatCounter {

  typedef typename Graph::EdgeId EdgeId;

 public:
  MatePairTransformStat(const Graph& g, const PairedInfoIndexT<Graph>& pair_info) :
        g_(g), pair_info_(pair_info), considered_dist_cnt_(0),
        unique_distance_cnt_(0), non_unique_distance_cnt_(0)
  {
  }

  virtual void Count() {
    for (auto it = pair_info_.begin(); it != pair_info_.end(); ++it) {
      de::Histogram infos = *it;
      EdgeId e1 = it.first();
      EdgeId e2 = it.second();
      ProcessInfos(e1, e2, infos);
    }
    INFO("Considered " << considered_dist_cnt_ << " edge pair distances (including trivial)");
    INFO(unique_distance_cnt_ << " edge distances connected with unique path of appropriate length");
    INFO(non_unique_distance_cnt_ << " edge distances connected with non-unique path of appropriate length");
  }

  size_t considered_edge_pair_count() {
    return considered_dist_cnt_;
  }

  size_t unique_distance_count() {
    return unique_distance_cnt_;
  }

  size_t non_unique_distance_count() {
    return non_unique_distance_cnt_;
  }

 private:
  const Graph& g_;
  const PairedInfoIndexT<Graph>& pair_info_;

  size_t considered_dist_cnt_;
  size_t unique_distance_cnt_;
  size_t non_unique_distance_cnt_;

  void ProcessInfos(EdgeId e1, EdgeId e2, const Histogram& infos) {
    for (auto it = infos.begin(); it != infos.end(); ++it) {
      Point point = *it;
      if (gr(point.d, 0.)) {
        if (eq(point.var, 0.)) {

          PathStorageCallback<Graph> counter(g_);

          PathProcessor<Graph> path_processor(g_,
              (size_t) (point.d - (double) g_.length(e1)),
              (size_t) (point.d - (double) g_.length(e1)),
              g_.EdgeEnd(e1), g_.EdgeStart(e2), counter);
          path_processor.Process();

          TRACE("Edges" << e1 << " : " << e2 << ": " << point.weight << " : " << point.d);
          TRACE("Path Numbs" << counter.size());

          if (counter.size() == 1)
            ++unique_distance_cnt_;
          if (counter.size() > 1)
            ++non_unique_distance_cnt_;
        }
        else
          non_unique_distance_cnt_++;

        considered_dist_cnt_++;
      }
    }
  }

  DECL_LOGGER("MatePairTransformStat")
};

template<class Graph>
class UniqueDistanceStat: public AbstractStatCounter {
  typedef omnigraph::de::PairedInfoIndexT<Graph> PairedIndex;

  const PairedIndex& paired_info_;
  size_t unique_;
  size_t non_unique_;
public:

  UniqueDistanceStat(const PairedIndex& paired_info) :
      paired_info_(paired_info), unique_(0), non_unique_(0) {

  }

  virtual ~UniqueDistanceStat() {

  }

  virtual void Count() {
    for (auto it = paired_info_.begin(); it != paired_info_.end(); ++it) {
      VERIFY((*it).size() > 0);
      if ((*it).size() > 1) {
        non_unique_++;
        //        for (auto info_it = (*it).begin(); info_it != (*it).end(); ++info_it) {
        //          //todo
        //        }
      } else {
        unique_++;
      }
    }INFO(unique_ << " unique edge distances");
    INFO(non_unique_ << " non unique edge distances");
  }

  size_t unique() {
    return unique_;
  }

  size_t non_unique() {
    return non_unique_;
  }
};

template<class Graph, class Index>
class EstimationQualityStat: public AbstractStatCounter {
private:
  typedef typename Graph::EdgeId EdgeId;
  typedef PairInfo<EdgeId> Info;
  typedef vector<Info> Infos;
  //input fields
  const Graph &graph_;
  const EdgeQuality<Graph, Index>& quality_;
  const PairedInfoIndex<Graph>& pair_info_;
  const PairedInfoIndex<Graph>& estimated_pair_info_;
  const PairedInfoIndex<Graph>& etalon_pair_info_;

  //output fields
  PairedInfoIndex<Graph> false_positives_;
  PairedInfoIndex<Graph> perfect_matches_;
  PairedInfoIndex<Graph> imperfect_matches_;
  PairedInfoIndex<Graph> false_negatives_;

//  PairedInfoIndexT<Graph> false_positive_weights_;
//  set<Info> false_positive_infos_;
//  vector<double> perfect_match_weights_;
////(weight, estimated_variance - actual_variance, number of etalon points)
//  vector<pair<pair<double, double> , size_t>> imperfect_match_stat_;
//  size_t false_negative_count_;
//  vector<Info> false_negative_infos_;

  bool CheckInterestInInfo(const Info& info) {
        if (math::ls(info.d, 0.)) return false;
        if (info.first == info.second && math::eq(info.d, 0.)) return false;
    return quality_.IsPositiveQuality(info.first)
            && quality_.IsPositiveQuality(info.second) && math::gr(info.weight, 0.);
  }

  void HandleFalsePositive(const Info& estimated) {
//    DEBUG("Handling false positive " << estimated);
    if (CheckInterestInInfo(estimated))
      false_positives_.AddPairInfo(estimated, false);
  }

  void HandleFalseNegative(const Info& etalon) {
    if (CheckInterestInInfo(etalon))
      false_negatives_.AddPairInfo(etalon, false);
  }

  void HandlePerfectMatch(const Info& etalon, const Info& estimated) {
    if (CheckInterestInInfo(estimated))
      perfect_matches_.AddPairInfo(estimated, false);
  }

  void HandleImperfectMatch(const Info &estimated_cluster,
      const Infos& etalon_matches) {
    if (CheckInterestInInfo(estimated_cluster))
      imperfect_matches_.AddPairInfo(estimated_cluster, false);
//    double etalon_variance = etalon_matches[etalon_matches.size() - 1].d
//        - etalon_matches[0].d;
//    imperfect_match_stat_.push_back(
//        make_pair(
//            make_pair(estimated_cluster.weight,
//                estimated_cluster.variance - etalon_variance),
//            etalon_matches.size()));
  }

//  void Flush() {
//    ProcessImperfectMatch(last_estimated_imperfect_match_,
//        last_etalon_imperfect_matches_);
//  }

  void HandlePairsNotInEtalon(
      const set<pair<EdgeId, EdgeId>>& pairs_in_etalon) {
    for (auto it = estimated_pair_info_.begin();
        it != estimated_pair_info_.end(); ++it) {
      Infos estimated_infos = *it;
      EdgeId first = estimated_infos[0].first;
      EdgeId second = estimated_infos[0].second;
      if (pairs_in_etalon.count(make_pair(first, second)) == 0) {
        //        for_each(estimated_infos.begin(), estimated_infos.end(),
        //            boost::bind(&EstimationQualityStat::HandleFalsePositive, this, _1));

        for (auto it2 = estimated_infos.begin();
            it2 != estimated_infos.end(); ++it2) {
          HandleFalsePositive(*it2);
        }
      }
    }
  }

  bool InfoLess(const Info& a, const Info& b) {
    if (eq(a.variance, 0.) && eq(b.variance, 0.)) {
      return ls(a.d, b.d);
    }
    return ls(a.d + a.variance, b.d - b.variance);
  }

  bool IsPerfectMatch(const Info& etalon, const Info& estimated) {
    return le(etalon.d, estimated.d) && ge(etalon.d, estimated.d)
        && eq(estimated.variance, 0.);
  }

  bool IsImperfectMatch(const Info& etalon, const Info& estimated) {
    return ge(etalon.d, estimated.d - estimated.variance)
        && le(etalon.d, estimated.d + estimated.variance);
  }

  size_t Move(size_t estimated_idx, const Infos &estimated_infos) {
    estimated_idx++;
    while (estimated_idx < estimated_infos.size()
        && math::eq(estimated_infos[estimated_idx].weight, 0.))
      estimated_idx++;
    return estimated_idx;
    return 0;
  }

  size_t InitIdx(const Infos &pair_infos) {
    return Move(-1, pair_infos);
  }

  void ProcessInfos(const Infos& etalon_infos, const Infos& estimated_infos) {
    size_t etalon_idx = InitIdx(etalon_infos);
    for (size_t estimated_idx = InitIdx(estimated_infos);
        estimated_idx < estimated_infos.size();
        estimated_idx = Move(estimated_idx, estimated_infos)) {
      while (estimated_idx < estimated_infos.size()
          && (etalon_idx == etalon_infos.size()
              || InfoLess(estimated_infos[estimated_idx],
                  etalon_infos[etalon_idx]))) {
        HandleFalsePositive(estimated_infos[estimated_idx]);
        estimated_idx = Move(estimated_idx, estimated_infos);
      }
      if (estimated_idx == estimated_infos.size()) {
        break;
      }
      while (etalon_idx < etalon_infos.size()
          && InfoLess(etalon_infos[etalon_idx],
              estimated_infos[estimated_idx])) {
        HandleFalseNegative(etalon_infos[etalon_idx]);
        etalon_idx = Move(etalon_idx, etalon_infos);
      }
      if (etalon_idx == etalon_infos.size()) {
        continue;
      }
      if (IsPerfectMatch(etalon_infos[etalon_idx],
          estimated_infos[estimated_idx])) {
        while (etalon_idx < etalon_infos.size()
            && IsPerfectMatch(etalon_infos[etalon_idx],
                estimated_infos[estimated_idx])) {
          HandlePerfectMatch(etalon_infos[etalon_idx],
              estimated_infos[estimated_idx]);
          etalon_idx = Move(etalon_idx, etalon_infos);
        }
      } else {
        vector<PairInfo<EdgeId> > cluster_hits;
        while (etalon_idx < etalon_infos.size()
            && IsImperfectMatch(etalon_infos[etalon_idx],
                estimated_infos[estimated_idx])) {
          cluster_hits.push_back(etalon_infos[etalon_idx]);
          etalon_idx = Move(etalon_idx, etalon_infos);
        }
        if (cluster_hits.size() == 0) {
          HandleFalsePositive(estimated_infos[estimated_idx]);
        } else {
          HandleImperfectMatch(estimated_infos[estimated_idx],
              cluster_hits);
        }
      }
    }
    //    for (size_t etalon_idx = 0; etalon_idx < etalon_infos.size(); ++etalon_idx) {
    //      Info etalon_info = etalon_infos[etalon_idx];
    ////      cout << "here" << endl;
    //      while (estimated_idx < estimated_infos.size() && InfoLess(estimated_infos[estimated_idx], etalon_info)) {
    //        HandleFalsePositive(estimated_infos[estimated_idx]);
    //        estimated_idx++;
    ////        cout << "here1" << endl;
    //      }
    ////      cout << "here2" << endl;
    //      if (estimated_idx != estimated_infos.size()
    //          && (HandleIfPerfectMatch(etalon_info, estimated_infos[estimated_idx])
    //              || HandleIfImperfectMatch(etalon_info, estimated_infos[estimated_idx]))) {
    //        last_matched = true;
    //      } else {
    //        HandleFalseNegative(etalon_info);
    //      }
    //    }
    //    if (last_matched)
    //      estimated_idx++;
    while (etalon_idx < etalon_infos.size()) {
      //      DEBUG("Handling false positives beyond all etalons");
      HandleFalseNegative(etalon_infos[etalon_idx]);
      etalon_idx = Move(etalon_idx, etalon_infos);
    }
    //    Flush();
  }

//  void ReportFalsePositiveWeights() {
//    sort(false_positive_weights_.begin(), false_positive_weights_.end());
//
//    INFO("False positive count: " << false_positive_weights_.size());
//  }
//
//  void ReportPerfectMatchWeights() {
//    sort(perfect_match_weights_.begin(), perfect_match_weights_.end());
//    INFO("Perfect match count: " << perfect_match_weights_.size());
//  }
//
//  void ReportImperfectMatchWeights() {
//    sort(imperfect_match_stat_.begin(), imperfect_match_stat_.end());
//    //todo do something better
//    INFO("Imperfect match count: " << imperfect_match_stat_.size());
//  }
//
//  void FalseNegativeCount() {
//    INFO("False negative count: " << false_negative_count_);
//  }

public:
  EstimationQualityStat(const Graph &graph,
      const EdgeQuality<Graph, Index>& quality,
      const PairedInfoIndex<Graph>& pair_info,
      const PairedInfoIndex<Graph>& estimated_pair_info,
      const PairedInfoIndex<Graph>& etalon_pair_info) :
      graph_(graph), quality_(quality), pair_info_(pair_info), estimated_pair_info_(
          estimated_pair_info), etalon_pair_info_(etalon_pair_info), false_positives_(
          graph_), perfect_matches_(graph_), imperfect_matches_(
          graph_), false_negatives_(graph_) {
  }

  virtual ~EstimationQualityStat() {
  }

  virtual void Count() {
    INFO("Counting distance estimation statistics");
    set<pair<EdgeId, EdgeId>> pairs_in_etalon;
    //    DEBUG("Handling pairs present in etalon information");
    for (auto it = etalon_pair_info_.begin(); it != etalon_pair_info_.end(); ++it) {
      Infos etalon_infos = *it;
      EdgeId first = etalon_infos[0].first;
      EdgeId second = etalon_infos[0].second;
      pairs_in_etalon.insert(make_pair(first, second));

      Infos estimated_infos = estimated_pair_info_.GetEdgePairInfo(first, second);
      //      DEBUG("Processing distances for pair " << first << ", " << second);
      ProcessInfos(etalon_infos, estimated_infos);
    }
    //    DEBUG("Handling pairs that are not in etalon information");
    HandlePairsNotInEtalon(pairs_in_etalon);

    INFO("FPR: " << fpr());
    INFO("FNR: " << fnr());
    INFO("Distance estimation statistics counted");
  }

  const PairedInfoIndexT<Graph>& false_positives() {
    return false_positives_;
  }

  const PairedInfoIndexT<Graph>& perfect_matches() {
    return perfect_matches_;
  }

  const PairedInfoIndexT<Graph>& imperfect_matches() {
    return imperfect_matches_;
  }

  const PairedInfoIndexT<Graph>& false_negatives() {
    return false_negatives_;
  }

  double fpr() {
    return 1. * false_positives_.size() / estimated_pair_info_.size();
  }

  double fnr() {
    return 1. * false_negatives_.size() / etalon_pair_info_.size();
  }

  void SaveStats(const string& dir_name) {
    //saving results
    INFO("Saving estimation statistic");
    make_dir(dir_name);
    typename PrinterTraits<Graph>::Printer printer(graph_);
    printer.savePaired(dir_name + "fp", false_positives_);
    printer.savePaired(dir_name + "pm", perfect_matches_);
    printer.savePaired(dir_name + "im", imperfect_matches_);
    printer.savePaired(dir_name + "fn", false_negatives_);
    INFO("Estimation statistics saved");
  }

//  vector<double> false_positive_weights() {
//    sort(false_positive_weights_.begin(), false_positive_weights_.end());
//    return false_positive_weights_;
//  }
//  vector<double> perfect_match_weights() {
//    sort(perfect_match_weights_.begin(), perfect_match_weights_.end());
//    return perfect_match_weights_;
//  }
//
//  vector<pair<pair<double, double> , size_t>> imperfect_match_weights() {
//    sort(imperfect_match_stat_.begin(), imperfect_match_stat_.end());
//    return imperfect_match_stat_;
//  }
//
//  size_t false_negative_count() {
//    return false_negative_count_;
//  }

//  void WriteFalseNegativeGaps(const string &file_name) {
//    ofstream stream;
//    stream.open(file_name);
//    vector<double> to_print;
//    //    for (size_t i = 0; i < false_negative_infos_.size(); i++) {
//    //      if (false_negative_infos_[i].d > 0)
//    //        to_print.push_back(
//    //            false_negative_infos_[i].d - graph_.length(
//    //                false_negative_infos_[i].first));
//    //    }
//    //    sort(to_print.begin(), to_print.end());
//    //    copy(to_print.begin(), to_print.end(),
//    //        ostream_iterator<double> (stream, "\n"));
//    for (size_t i = 0; i < false_negative_infos_.size(); i++) {
//      stream << false_negative_infos_[i] << endl;
//    }
//    stream.close();
//  }
//
//  void WriteEstmationStats(const string &output_folder) {
//    ofstream stream;
//    stream.open(output_folder + "/perfect.inf");
//    copy(perfect_match_weights_.begin(), perfect_match_weights_.end(),
//        ostream_iterator<double>(stream, "\n"));
//    stream.close();
//
//    stream.open(output_folder + "/false_positive.inf");
//    copy(false_positive_weights_.begin(), false_positive_weights_.end(),
//        ostream_iterator<double>(stream, "\n"));
//    stream.close();
//    WriteWorstEdgesStat(output_folder, 1000000);
//  }

//  void WriteEdgePairInfo(const string &file_name, Infos infos) {
//    ofstream stream;
//    stream.open(file_name);
//    for (size_t i = 0; i < infos.size(); i++) {
//      stream << infos[i] << endl;
//    }
//    stream.close();
//  }
//
//  string ConstructEdgePairFileName(const string output_folder,
//      const string &name, const string &modifier, size_t index) {
//    stringstream ss;
//    ss.clear();
//    ss << output_folder << "/" << name << "_" << index << "_" << modifier
//        << ".inf";
//    return ss.str();
//  }

//  void WriteWorstEdgesStat(const string &output_folder, double bound) {
//    size_t count = 0;
//    WriteFalseNegativeGaps(output_folder + "/gaps.inf");
//    for (auto iterator = false_positive_infos_.begin();
//        iterator != false_positive_infos_.end(); ++iterator) {
//      if (iterator->weight > bound) {
//        WriteEdgePairInfo(
//            ConstructEdgePairFileName(output_folder, "fp",
//                "histogram", count),
//            pair_info_.GetEdgePairInfo(iterator->first,
//                iterator->second));
//        WriteEdgePairInfo(
//            ConstructEdgePairFileName(output_folder, "fp",
//                "estimated", count),
//            estimated_pair_info_.GetEdgePairInfo(iterator->first,
//                iterator->second));
//        WriteEdgePairInfo(
//            ConstructEdgePairFileName(output_folder, "fp", "etalon",
//                count),
//            etalon_pair_info_.GetEdgePairInfo(iterator->first,
//                iterator->second));
//        count++;
//      }
//    }
//    for (auto iterator = false_negative_infos_.begin();
//        iterator != false_negative_infos_.end(); ++iterator) {
//      if (iterator->weight > bound) {
//        WriteEdgePairInfo(
//            ConstructEdgePairFileName(output_folder, "fp",
//                "histogram", count),
//            pair_info_.GetEdgePairInfo(iterator->first,
//                iterator->second));
//        WriteEdgePairInfo(
//            ConstructEdgePairFileName(output_folder, "fp",
//                "estimated", count),
//            estimated_pair_info_.GetEdgePairInfo(iterator->first,
//                iterator->second));
//        WriteEdgePairInfo(
//            ConstructEdgePairFileName(output_folder, "fp", "etalon",
//                count),
//            etalon_pair_info_.GetEdgePairInfo(iterator->first,
//                iterator->second));
//        count++;
//      }
//    }
//  }

};

template<class Graph>
class ClusterStat: public AbstractStatCounter {

  typedef typename Graph::EdgeId EdgeId;
  typedef pair<double, double> DoublePair;

 public:
  ClusterStat(const PairedInfoIndexT<Graph>& estimated_pair_info) :
    estimated_pair_info_(estimated_pair_info)
  {
  }

  virtual ~ClusterStat()
  {
  }

  virtual void Count() {
    for (auto it = estimated_pair_info_.begin(); it != estimated_pair_info_.end(); ++it) {
      de::Histogram infos = *it;
      for (auto it2 = infos.begin(); it2 != infos.end(); ++it2) {
        Point point = *it2;
        if (gr(point.var, 0.))
          weight_variance_stat_.push_back(make_pair(point.weight, point.var));
      }
      //todo talk with Anton!!!
      //      for (auto it2 = (*it).begin(); it2 != (*it).end(); ++it2) {
      //        Info info = *it2;
      ////        if (gr(info.variance, 0)) {
      //          weight_variance_stat_.push_back(make_pair(info.weight, info.variance));
      ////        }
      //      }
    }
    stringstream ss;
    copy(weight_variance_stat_.begin(), weight_variance_stat_.end(),
        ostream_iterator<DoublePair>(ss, ", "));
    INFO("Estimated cluster stat: " << ss.str());
  }

  vector<DoublePair> weight_variance_stat() {
    sort(weight_variance_stat_.begin(), weight_variance_stat_.end());
    return weight_variance_stat_;
  }

private:
  const PairedInfoIndexT<Graph>& estimated_pair_info_;
  vector<DoublePair> weight_variance_stat_;

  DECL_LOGGER("EstimatedClusterStat");
};

template<class Graph, class Index>
class CoverageStatistics{

private:
    Graph& graph_;
    EdgeQuality<Graph, Index> & edge_qual_;

    bool DoesSuit(VertexId vertex){
        bool ans = true;
        for (size_t i = 0; ans && i<graph_.OutgoingEdgeCount(vertex); i++)
            ans = ans & math::gr(edge_qual_.quality(graph_.OutgoingEdges(vertex)[i]), 0.);
        for (size_t i = 0; ans && i<graph_.IncomingEdgeCount(vertex); i++)
            ans = ans & math::gr(edge_qual_.quality(graph_.IncomingEdges(vertex)[i]), 0.);
        return ans;
    }

public:
    CoverageStatistics(Graph& graph, EdgeQuality<Graph, Index>& edge_qual):
    graph_(graph), edge_qual_(edge_qual){
    }

    virtual ~CoverageStatistics(){}

    virtual void Count(){

        map<double, size_t> cov_map;
        map<double, size_t> ratio_map;
        map<double, size_t> len_map;
        size_t area = 0;
        size_t area15 = 0;
        size_t area10 = 0;
        size_t area5 = 0;
        size_t area2 = 0;
        for (auto iter = graph_.ConstEdgeBegin(); !iter.IsEnd(); ++iter){
            len_map[graph_.length(*iter)]++;
        }
        for (auto iter = graph_.begin(); iter != graph_.end(); ++iter)
            if (true || DoesSuit(*iter) ){

                double plus_cov = 0.;
                double min_cov = 0.;
                double plus_all_cov = 0.;
                double min_all_cov = 0.;
                bool suit_us = true;

                if (graph_.IncomingEdgeCount(*iter)*graph_.OutgoingEdgeCount(*iter) == 0) continue;

                for (size_t i = 0; suit_us && i<graph_.IncomingEdgeCount(*iter); i++)
                    if (graph_.length(graph_.IncomingEdges(*iter)[i]) < 80){
                        if (math::ge(edge_qual_.quality(graph_.IncomingEdges(*iter)[i]), 1.))
                            plus_cov += graph_.coverage(graph_.IncomingEdges(*iter)[i]);
                        plus_all_cov += graph_.coverage(graph_.IncomingEdges(*iter)[i]);
                    }else suit_us = false;
                for (size_t i = 0; suit_us && i<graph_.OutgoingEdgeCount(*iter); i++)
                    if (graph_.length(graph_.OutgoingEdges(*iter)[i]) < 80){
                        if (math::ge(edge_qual_.quality(graph_.OutgoingEdges(*iter)[i]), 1.))
                            min_cov += graph_.coverage(graph_.OutgoingEdges(*iter)[i]);
                        min_all_cov += graph_.coverage(graph_.OutgoingEdges(*iter)[i]);
                    }else suit_us = false;

                if (!suit_us) continue;

                if (math::eq(min_cov, 0.) || math::eq(plus_cov, 0.)) continue;

                double delta_cov = math::round(1000.*(plus_cov - min_cov)/(plus_cov + min_cov));

                double ratio_cov = math::round(1000.*(plus_cov + min_cov)/(plus_all_cov + min_all_cov));

                if (math::ls(abs(delta_cov), 150.)) area15++;
                if (math::ls(abs(delta_cov), 100.)) area10++;
                if (math::ls(abs(delta_cov), 50.)) area5++;
                if (math::ls(abs(delta_cov), 20.)) area2++;
                area++;

                cov_map[delta_cov/10.]++;
                ratio_map[ratio_cov/10.]++;

        }

        for (auto iter = ratio_map.begin(); iter != ratio_map.end(); ++iter){
            INFO("Ratio " << (*iter).first << " " << (*iter).second);
        }

        for (auto iter = cov_map.begin(); iter != cov_map.end(); ++iter){
            INFO("Cov " << (*iter).first << " " << (*iter).second);
        }

        INFO("stats_cov "  << area << " " << area2 << " " << area5 << " " << area10 << " " << area15);

        for (auto iter = len_map.begin(); iter != len_map.end(); ++iter){
            INFO("Len " << (*iter).first << " " << (*iter).second);
        }

    }

};
