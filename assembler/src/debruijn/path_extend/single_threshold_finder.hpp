/*
 * single_threshold_finder.hpp
 *
 *  Created on: Mar 1, 2013
 *      Author: ira
 */

#ifndef SINGLE_THRESHOLD_FINDER_HPP_
#define SINGLE_THRESHOLD_FINDER_HPP_
#include "path_extend/path_extend_launch.hpp"
#include "path_extend/paired_library.hpp"
#include "late_pair_info_count.hpp"
#include "graphio.hpp"
using namespace debruijn_graph;
class SingleThresholdFinder {

public:
	int insert_size_min;
	int insert_size_max;
	int length_to_split;
	bool is_mp_;
	SingleThresholdFinder(int insert_size_min1, int insert_size_max1, int length_to_split1, bool is_mp = false):insert_size_min(insert_size_min1),
			insert_size_max(insert_size_max1), length_to_split(length_to_split1), is_mp_(is_mp) {
	}

	double find_threshold() {
		Sequence genome = cfg::get().developer_mode ? cfg::get().ds.reference_genome : Sequence();
		conj_graph_pack gp(cfg::get().K, cfg::get().output_dir, genome,
				cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling,
				!cfg::get().developer_mode);
		INFO("single threshold finder");
		exec_simplification(gp);
		gp.index.Detach();
		set<BidirectionalPath*> goodPaths;
		split_long_edges(gp, length_to_split, goodPaths);
		gp.index.Attach();
		gp.index.Refill();
		PairedIndexT paired_index(gp.g);
		size_t is = cfg::get().ds.IS();
		io::ReadStreamVector<io::IReader<io::PairedReadSeq>> paired_streams = paired_binary_readers(true, is);
		FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index,
				gp.kmer_mapper, paired_index, paired_streams, gp.k_value);
		PairedIndexT clustered_index(gp.g);
		if (!is_mp_){
			estimate_distance(gp, paired_index, clustered_index);
		}
		PairedInfoLibrary* lib_not_cl = new PairedInfoLibrary(cfg::get().K, gp.g, cfg::get().ds.RL(), is, cfg::get().ds.is_var(), paired_index);
		PairedInfoLibrary* lib_cl = new PairedInfoLibrary(cfg::get().K, gp.g, cfg::get().ds.RL(), is, cfg::get().ds.is_var(), clustered_index);

		map<PairInfo<EdgeId>, double> good_pi;
		map<PairInfo<EdgeId>, double> bad_pi;
		for (auto iter = goodPaths.begin(); iter != goodPaths.end(); ++iter) {
			analyze_one_path(gp, *iter, lib_not_cl, lib_cl, good_pi, bad_pi);
		}
		writeToFile(gp, good_pi, bad_pi, lib_not_cl);
		deletePaths(goodPaths);
		vector<double> good_pi_val;
		get_norm_pi(gp, good_pi, good_pi_val);
		vector<double> bad_pi_val;
		get_norm_pi(gp, bad_pi, bad_pi_val);
		return find_intersection(good_pi_val, bad_pi_val);
	}

private:
	string edge_info(conj_graph_pack& gp, string pref,
			const PairInfo<EdgeId>& pi, double ideal_pi,
			PairedInfoLibrary* lib) {
		stringstream ss;
		double w = pi.point.weight;
		double d = pi.point.d;
		double cov1 = gp.g.coverage(pi.first);
		double cov2 = gp.g.coverage(pi.second);
		double pi1 = lib->get_all_pi_count(pi.first, true) / gp.g.length(pi.first);
		double pi2 = lib->get_all_pi_count(pi.second, false) / gp.g.length(pi.second);
		double pi_norm1 = lib->get_all_norm_pi(pi.first, true);
		double pi_norm2 = lib->get_all_norm_pi(pi.second, false);
		double pi_norm1_aver = lib->get_all_norm_pi_aver(pi.first, true);
		double pi_norm2_aver = lib->get_all_norm_pi_aver(pi.second, false);
		double test1 = w / ideal_pi;
		double test2 = (w / ideal_pi) / min(pi1, pi2);
		double test3 = (w / ideal_pi) / min(pi_norm1, pi_norm2);
		double test4 = test1 / min(pi_norm1_aver, pi_norm2_aver);
		ss << pref << " " << (w / ideal_pi) << " " << w << " " << ideal_pi
				<< " " << gp.g.str(pi.first) << " " << gp.g.str(pi.second)<<
				" " << test4
				<< " " << d << " " << cov1 << " " << cov2 << " "
				<< (w / ideal_pi) / min(cov1, cov2) << " " << pi1 << " "
				<< (w / ideal_pi) / pi1 << " " << (w / ideal_pi) / min(pi1, pi2)
				<< " " << (w / ideal_pi) / min(pi_norm1, pi_norm2) << "\n";
		return ss.str();
	}

	void writeToFile(conj_graph_pack& gp, map<PairInfo<EdgeId>, double>& good_pi,  map<PairInfo<EdgeId>, double>& bad_pi, PairedInfoLibrary* lib){
		string file_name = cfg::get().output_dir + "single_threshold.txt";
		ofstream os(file_name);
		os << "good pi " << good_pi.size() << " bad pi "<<bad_pi.size() << "\n";
		for (auto iter = bad_pi.begin(); iter != bad_pi.end(); ++iter) {
			os << edge_info(gp, " all_bad", iter->first, iter->second, lib);
		}
		for (auto iter = good_pi.begin(); iter != good_pi.end(); ++iter) {
			os << edge_info(gp, " good", iter->first, iter->second, lib);
		}
		os.close();
	}

	void get_norm_pi(conj_graph_pack& gp, map<PairInfo<EdgeId>, double>& pi, vector<double>& result){
		for (auto iter = pi.begin(); iter != pi.end(); ++iter) {
			if (iter->first.point.weight > 0.0 and iter->second > 0.0 and
					!boost::math::isnan<double>(iter->first.point.weight / iter->second) and
					!boost::math::isinf<double>(iter->first.point.weight / iter->second)){
				result.push_back(iter->first.point.weight / iter->second);
			}
		}
	}

	double find_intersection(vector<double>& good_pi, vector<double>& bad_pi) {
		std::sort(good_pi.begin(), good_pi.end());
		std::sort(bad_pi.begin(), bad_pi.end());
		size_t good_iter = 0;
		size_t bad_iter = 0;
		double cur_threshold = 0.0;
		double good_percent = 0;
		double bad_percent = 1;
		while (good_percent < bad_percent and good_iter < good_pi.size() and bad_iter < bad_pi.size()) {
			cur_threshold = good_pi[good_iter];
			while (bad_iter < bad_pi.size() and bad_pi[bad_iter] <= cur_threshold) {
				bad_iter++;
			}
			good_percent = (double)good_iter / (double)good_pi.size();
			bad_percent = 1 - (double)bad_iter / (double)bad_pi.size();
			good_iter += 1;
		}
		return cur_threshold;
	}


	void deletePaths(set<BidirectionalPath*>& paths) {
		for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
			delete (*iter);
		}
	}

	void split_long_edges(conj_graph_pack& gp, size_t length, set<BidirectionalPath*>& result) {
		INFO("begin split long edges");
		set<EdgeId> longEdges;
		Graph& graph = gp.g;
		size_t com_length = 0;
		for (auto edge = graph.SmartEdgeBegin(); !edge.IsEnd(); ++edge) {
			if (graph.length(*edge) > length && longEdges.find(gp.g.conjugate(*edge)) == longEdges.end()) {
				longEdges.insert(*edge);
				com_length += graph.length(*edge);
			}
		}
		for (auto longEdge = longEdges.begin(); longEdge != longEdges.end(); ++longEdge){
			if (graph.int_id(*longEdge) <= 0){
				continue;
			}
			vector<EdgeId> edge_path;
			EdgeId secondEdge = *longEdge;
			while (graph.length(secondEdge) > length) {
				pair<EdgeId, EdgeId> result_edges = graph.SplitEdge(secondEdge, length);
				edge_path.push_back(result_edges.first);
				secondEdge = result_edges.second;
			}
			edge_path.push_back(secondEdge);
			result.insert(new BidirectionalPath(graph, edge_path));
		}
		INFO("end split long edges");
	}

	void find_idel_pair_info(conj_graph_pack& gp, size_t from, BidirectionalPath* path,  PairedInfoLibrary* lib, map<EdgeId, PairInfo<EdgeId> >& idealPairInfo){
		size_t length = 0;
		for (size_t index = from; index < path->Size() && length < (size_t)insert_size_max; ++index){
			double ideal_pi_w = lib->IdealPairedInfo(path->At(from), path->At(index), length);
			if (ideal_pi_w > 0){
				PairInfo<EdgeId> ideal_pi(path->At(from), path->At(index), length, ideal_pi_w, 0);
				idealPairInfo.insert(make_pair(path->At(index), ideal_pi));
			}
			length += gp.g.length(path->At(index));
		}
	}

	void analyze_one_path(conj_graph_pack& gp, BidirectionalPath* path, PairedInfoLibrary* lib_not_cl,PairedInfoLibrary* lib_cl,
			map<PairInfo<EdgeId>, double >& good_pi, map<PairInfo<EdgeId>, double>& bad_pi){
		size_t length_begin_path = 0;
		for (size_t index = 0; index < path->Size() && length_begin_path < path->Length() - insert_size_max; ++index){
			length_begin_path += gp.g.length(path->At(index));

			vector<PairInfo<EdgeId> > edge_pi= lib_not_cl->index_.GetEdgeInfo(path->At(index));
			vector<PairInfo<EdgeId> > edge_pi_cl= lib_cl->index_.GetEdgeInfo(path->At(index));
			map<EdgeId, PairInfo<EdgeId> > idealPairInfo;
			set<EdgeId> good_edges;
			find_idel_pair_info(gp, index, path, lib_not_cl, idealPairInfo);
			if (!is_mp_){
				for (size_t i_pi = 0; i_pi < edge_pi_cl.size(); ++i_pi){
					if (idealPairInfo.find(edge_pi_cl[i_pi].second) != idealPairInfo.end()){
						good_pi[edge_pi_cl[i_pi]] = idealPairInfo.find(edge_pi_cl[i_pi].second)->second.point.weight;
						good_edges.insert(edge_pi_cl[i_pi].second);
					}
				}
			}
			map<pair<EdgeId, EdgeId>, vector<Point> > good_edge_pi;
			map<pair<EdgeId, EdgeId>, vector<Point> > bad_edge_pi;
			for (size_t i_pi = 0; i_pi < edge_pi.size(); ++i_pi){
				if (idealPairInfo.find(edge_pi[i_pi].second)== idealPairInfo.end()){
					add_pi(edge_pi[i_pi], bad_edge_pi);
				} else {
					add_pi(edge_pi[i_pi], good_edge_pi);
				}
			}
			if (is_mp_){
				for (auto iter = good_edge_pi.begin(); iter != good_edge_pi.end(); ++iter){
					double weight = lib_not_cl->find_weight(iter->second);
					if (weight > 0.0){
						PairInfo<EdgeId> ideal_pi = idealPairInfo.find(iter->first.second)->second;
						PairInfo<EdgeId> real_pi(ideal_pi.first, ideal_pi.second, ideal_pi.point.d, weight, ideal_pi.point.var);
						good_pi[real_pi] = ideal_pi.point.weight;
					}
				}
			}
			lib_not_cl->clust_pi(bad_edge_pi, bad_pi);
		}
	}

	void add_pi(PairInfo<EdgeId>& pi, map<pair<EdgeId, EdgeId>, vector<Point> >& all_pi) {
		if (pi.point.d > 0) {
			auto pair_info = make_pair(pi.first, pi.second);
			if (all_pi.find(pair_info) == all_pi.end()) {
				all_pi.insert(make_pair(pair_info, vector<Point>()));
			}
			all_pi[pair_info].push_back(pi.point);
		}
	}
};


#endif /* SINGLE_THRESHOLD_FINDER_HPP_ */
