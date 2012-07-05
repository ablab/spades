/*
 * pair_info_improver.hpp
 *
 *  Created on: Jul 4, 2012
 *      Author: avsirotkin
 */

#pragma once

#include "standard.hpp"
#include "path_utils.hpp"
#include "graph_pack.hpp"
#include "split_path_constructor.hpp"



namespace debruijn_graph{


template<class Graph>
int DeleteIfExist(const Graph& g, PairedInfoIndex<Graph>& clustered_index,
		const PairInfo<typename Graph::EdgeId>& p_info) {
	int cnt = 0;
	auto pi_vector = clustered_index.GetEdgePairInfo(p_info.first, p_info.second);
	for (auto pi_iter = pi_vector.begin(); pi_iter != pi_vector.end();
			++pi_iter) {
		if (abs(pi_iter->d - p_info.d) < 0.01) {
			clustered_index.RemovePairInfo(*pi_iter);
			clustered_index.RemovePairInfo(BackwardInfo(*pi_iter));
			cnt += 2;
		}
	}
	TRACE("cnt += "<<cnt);
	return cnt;
}



template<class Graph>
int DeleteConjugatePairInfo(const Graph& g, PairedInfoIndex<Graph>& clustered_index,
		const PairInfo<typename Graph::EdgeId>& p_info) {
	int cnt = 0;
	auto pi_vector = clustered_index.GetEdgePairInfo(
			g.conjugate(p_info.second),
			g.conjugate(p_info.first));
	double tmpd = p_info.d + g.length(p_info.second)
			- g.length(p_info.first);
	for (auto pi_iter = pi_vector.begin(); pi_iter != pi_vector.end();
			++pi_iter) {
//					INFO("dist "<< tmpd << " versus "<<pi_iter->d<<" var "<<pi_iter->variance);
		if (abs(pi_iter->d - tmpd) < 0.01) {
			clustered_index.RemovePairInfo(*pi_iter);
			clustered_index.RemovePairInfo(BackwardInfo(*pi_iter));
			cnt += 2;
			DEBUG("Removed pi "<< g.int_id(pi_iter->first) << " "<< g.int_id(pi_iter->second)<<" d "<<pi_iter->d<<" var "<<pi_iter->variance);
		}
	}
	TRACE("cnt += "<<cnt);
	return cnt;
}




template<class Graph>
bool TryToAddPairInfo(const Graph& g,
		PairedInfoIndex<Graph>& clustered_index,
		const PairInfo<typename Graph::EdgeId>& p_i) {
		return TryToAddPairInfo(g, clustered_index, p_i.first, p_i.second, p_i.d, p_i.weight, p_i.variance);
}



template<class Graph>
bool TryToAddPairInfo(const Graph &g,
		PairedInfoIndex<Graph>& clustered_index,
		typename Graph::EdgeId first,
		typename Graph::EdgeId second, double tmpd, double w, double var = 0) {
	auto pi_vector = clustered_index.GetEdgePairInfo(first, second);
	bool already_exist = false;
	for (auto pi_iter = pi_vector.begin(); pi_iter != pi_vector.end();
			++pi_iter) {
		if (abs(pi_iter->d - tmpd) < 0.01 + pi_iter->variance + var) {
			already_exist = true;
		}
	}
	if (!already_exist) {
		clustered_index.AddPairInfo(
				PairInfo<typename Graph::EdgeId>(first, second,
						tmpd, w, var));
		clustered_index.AddPairInfo(
				PairInfo<typename Graph::EdgeId>(
						g.conjugate(second),
						g.conjugate(first),
						tmpd - g.length(first) + g.length(second), w, var));
		return true;
	}
	return false;
}





template<class Graph>
class PairInfoInprover{
	const Graph &g_;
public:
	PairInfoInprover(const Graph& g): g_(g) {};
	void ImprovePairedInfo(PairedInfoIndex<Graph>& clustered_index, bool parallel = false, size_t num_treads = 1){
		if (parallel) {
			ParalelCorrectPairedInfo(clustered_index, num_treads);
			ParalelCorrectPairedInfo(clustered_index, num_treads);
		} else {
			NonParalelCorrectPairedInfo(clustered_index);
			NonParalelCorrectPairedInfo(clustered_index);

		}
	}

private:
	void ParalelCorrectPairedInfo(PairedInfoIndex<Graph>& clustered_index, const size_t nthreads = 4) {
		int missing_paired_info_count = 0;
		int extra_paired_info_count = 0;

		extra_paired_info_count = ParalelRemoveContraditional(clustered_index, nthreads);
		missing_paired_info_count = ParalelFillMissing(clustered_index, nthreads);

		INFO(
				"Paired info stats: missing = " << missing_paired_info_count
						<< "; contradictional = " << extra_paired_info_count);
	}

	void NonParalelCorrectPairedInfo(PairedInfoIndex<Graph>& clustered_index, const size_t nthreads = 4) {
		int missing_paired_info_count = 0;
		int extra_paired_info_count = 0;

		extra_paired_info_count = NonParalelRemoveContraditional(clustered_index);
		missing_paired_info_count = NonParalelFillMissing(clustered_index);

		INFO(
				"Paired info stats: missing = " << missing_paired_info_count
						<< "; contradictional = " << extra_paired_info_count);
	}



	int ParalelRemoveContraditional(PairedInfoIndex<Graph>& clustered_index, size_t nthreads = 4){

		int cnt = 0;
		DEBUG("ParalelRemoveContraditional: Put infos to vector");
		std::vector<std::vector<PairInfo<typename Graph::EdgeId>>> infos;
		for (auto e_iter = g_.SmartEdgeBegin(); !e_iter.IsEnd();
				++e_iter) {
			if (g_.length(*e_iter)	>= cfg::get().rr.max_repeat_length)
				infos.push_back(clustered_index.GetEdgeInfo(*e_iter));
		}

		std::vector< PairInfoIndexData<typename Graph::EdgeId>* > to_remove(nthreads);
	    for (size_t i = 0; i < nthreads; ++i) {
	        to_remove[i] = new PairInfoIndexData<typename Graph::EdgeId>();
	    }
	    DEBUG("ParalelRemoveContraditional: Start threads");
	//	size_t n = 0;
		#pragma omp parallel num_threads(nthreads)
		{
			#pragma omp for
			for (size_t i = 0; i < infos.size(); ++i)
			{
				FindInconsistent(clustered_index, infos[i], to_remove[omp_get_thread_num()]);
	//            VERBOSE_POWER_T(++n, 100, " edge pairs processed. Cur thread "<<omp_get_thread_num());

			}
		}

		DEBUG("ParalelRemoveContraditional: Threads finished");

		for (size_t i = 0; i < nthreads; ++i) {
			for (auto remove_iter = (*to_remove[i]).begin(); remove_iter != (*to_remove[i]).end(); ++remove_iter) {
				cnt += DeleteIfExist<Graph>(g_, clustered_index, *remove_iter);
				cnt += DeleteConjugatePairInfo<Graph>(g_, clustered_index, *remove_iter);
			}
	        delete to_remove[i];
		}
		DEBUG("ParalelRemoveContraditional: Clean finished");
		return cnt;
	}

	int NonParalelRemoveContraditional(PairedInfoIndex<Graph>& clustered_index){

		int cnt = 0;

		PairInfoIndexData<typename Graph::EdgeId> *to_remove = new PairInfoIndexData<typename Graph::EdgeId>();

		for (auto e_iter = g_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
			std::vector<PairInfo<typename Graph::EdgeId>> edge_infos = clustered_index.GetEdgeInfo(*e_iter);
			if (g_.length(*e_iter)	>= cfg::get().rr.max_repeat_length)
				FindInconsistent(clustered_index, edge_infos, to_remove);
		}

		for (auto remove_iter = (*to_remove).begin(); remove_iter != (*to_remove).end(); ++remove_iter) {
			cnt += DeleteIfExist<Graph>(g_, clustered_index, *remove_iter);
			cnt += DeleteConjugatePairInfo<Graph>(g_, clustered_index, *remove_iter);
		}

		delete to_remove;
		return cnt;
	}


	void FindInconsistent(PairedInfoIndex<Graph>& clustered_index, std::vector<PairInfo<typename Graph::EdgeId>>& edge_infos, PairInfoIndexData<typename Graph::EdgeId>* pi) {

		for(size_t i = 0; i < edge_infos.size(); i++) {
			for(size_t j = 0; j < edge_infos.size(); j++) {
				if (i!=j){
					if (!isConsistent(edge_infos[i], edge_infos[j])){
						TRACE("Inconsistent!!!")
						if (edge_infos[i].weight > edge_infos[j].weight) {
							pi->AddPairInfo(edge_infos[j]);
						} else  {
							pi->AddPairInfo(edge_infos[i]);
						}
					}
				}
			}
		}
	}


	bool isConsistent(PairInfo<typename Graph::EdgeId>& first_info,
			PairInfo<typename Graph::EdgeId>& second_info) {
		if ((first_info.d < 0.0001 || second_info.d < 0.0001)
						|| (first_info.d > second_info.d)) return true;
	//	size_t max_comparable_path = *cfg::get().ds.IS - K
	//			+ size_t(*cfg::get().ds.is_var);
		int pi_distance = second_info.d - first_info.d;
		int first_length = g_.length(first_info.second);
		double variance = first_info.variance + second_info.variance;

		auto first_edge = first_info.second;
		auto second_edge = second_info.second;
		TRACE("   PI "<< g_.int_id(first_info.first)<<" "<<g_.int_id(first_info.second)<<" d "<<first_info.d<<"var "<<first_info.variance<<" tr "<< omp_get_thread_num());
		TRACE("vs PI "<< g_.int_id(second_info.first)<<" "<<g_.int_id(second_info.second)<<" d "<<second_info.d<<"var "<<second_info.variance<<" tr "<< omp_get_thread_num());


		if (abs(pi_distance - first_length)<=variance){
			if (g_.EdgeEnd(first_edge) == g_.EdgeStart(second_edge)) return true;
			else {
				auto paths = GetAllPathsBetweenEdges(g_, first_edge, second_edge, 0,
						pi_distance - first_length + variance);
				return (paths.size() > 0);
			}
		}
		else {
			if ((int)pi_distance > first_length){
				auto paths = GetAllPathsBetweenEdges(g_, first_edge, second_edge, pi_distance - first_length - variance,
						pi_distance - first_length + variance);
				return (paths.size() > 0);
			}
			return false;
		}
	}



	int ParalelFillMissing(PairedInfoIndex<Graph>& clustered_index, size_t nthreads = 4){
		int cnt = 0;
		DEBUG("Fill missing: Put infos to vector");
		std::vector<std::vector<PairInfo<typename Graph::EdgeId>>> infos;
		for (auto e_iter = g_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
				infos.push_back(clustered_index.GetEdgeInfo(*e_iter));
		}

		std::vector< PairedInfoIndex<Graph>* > to_add(nthreads);
	    for (size_t i = 0; i < nthreads; ++i) {
	        to_add[i] = new PairedInfoIndex<Graph>(g_);
	    }
	    SplitPathConstructor<Graph> spc(g_);
	    DEBUG("Fill missing: Start threads");
	//	size_t n = 0;
		#pragma omp parallel num_threads(nthreads)
		{
			#pragma omp for
			for (size_t i = 0; i < infos.size(); ++i)
			{
	//			FindInconsistent<graph_pack>(origin_gp, clustered_index, infos[i], to_add[omp_get_thread_num()]);
				vector<PathInfoClass<Graph>> paths = spc.ConvertEdgePairInfoToSplitPathes(infos[i]);
				for (auto iter = paths.begin(); iter != paths.end(); iter ++) {
					DEBUG("Path "<<iter->PrintPath(g_));
					for (auto pi_iter = iter->begin(); pi_iter != iter->end(); pi_iter ++) {
						TryToAddPairInfo(g_, *to_add[omp_get_thread_num()], *pi_iter);
					}
				}
	//			VERBOSE_POWER_T(++n, 100, " edge pairs processed. Cur thread "<<omp_get_thread_num());

			}
		}

		DEBUG("Fill missing: Threads finished");

		for (size_t i = 0; i < nthreads; ++i) {
			for (auto add_iter = (*to_add[i]).begin(); add_iter != (*to_add[i]).end(); ++add_iter) {
				auto pinfos = *add_iter;
				for (auto add_pi_iter = pinfos.begin(); add_pi_iter != pinfos.end(); ++add_pi_iter) {
					if (TryToAddPairInfo(g_, clustered_index, *add_pi_iter))
					{
						DEBUG("Fill missing: PI added "<<cnt);
						cnt++;
					}
				}
			}
	        delete to_add[i];
		}
		DEBUG("Fill missing: Clean finished");
		return cnt;
	}


	int NonParalelFillMissing(PairedInfoIndex<Graph>& clustered_index){
		int cnt = 0;


		PairedInfoIndex<Graph> to_add(g_);
	    SplitPathConstructor<Graph> spc(g_);
		for (auto e_iter = g_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
			std::vector<PairInfo<typename Graph::EdgeId>> edge_infos = clustered_index.GetEdgeInfo(*e_iter);
			vector<PathInfoClass<Graph>> paths = spc.ConvertEdgePairInfoToSplitPathes(edge_infos);
			for (auto iter = paths.begin(); iter != paths.end(); iter ++) {
				DEBUG("Path "<<iter->PrintPath(g_));
				for (auto pi_iter = iter->begin(); pi_iter != iter->end(); pi_iter ++) {
					TryToAddPairInfo(g_, to_add, *pi_iter);
				}
			}
		}


		for (auto add_iter = (to_add).begin(); add_iter != (to_add).end(); ++add_iter) {
			auto pinfos = *add_iter;
			for (auto add_pi_iter = pinfos.begin(); add_pi_iter != pinfos.end(); ++add_pi_iter) {
				if (TryToAddPairInfo(g_, clustered_index, *add_pi_iter))
				{
					DEBUG("Fill missing: PI added "<<cnt);
					cnt++;
				}
			}
		}
		return cnt;
	}





	DECL_LOGGER("PairInfoImprover")
};



//template<class graph_pack>
//void CorrectPairedInfo(const graph_pack& origin_gp,
//		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index,
//		bool clean = true, bool add = true) {
//	size_t k = graph_pack::k_value;
//	size_t delta = size_t(*cfg::get().ds.is_var);
//	size_t max_comparable_path = *cfg::get().ds.IS + delta - k;
//	int missing_paired_info_count = 0;
//	int extra_paired_info_count = 0;
//	int long_edges_count = 0;
//
//	DEBUG("Using max path cutoff = " << max_comparable_path);
//	for (auto e_iter = origin_gp.g.SmartEdgeBegin(); !e_iter.IsEnd();
//			++e_iter) {
//		if (origin_gp.g.length(*e_iter) >= cfg::get().rr.max_repeat_length)
//			long_edges_count++;
//		auto pi = clustered_index.GetEdgeInfo(*e_iter);
//		for (auto i_iter = pi.begin(); i_iter != pi.end(); ++i_iter) {
//			for (auto j_iter = pi.begin(); j_iter != pi.end(); ++j_iter) {
//				if (i_iter != j_iter) {
//					PairInfo<typename graph_pack::graph_t::EdgeId> first_info =
//							*i_iter;
//					PairInfo<typename graph_pack::graph_t::EdgeId> second_info =
//							*j_iter;
//					if (origin_gp.g.length(*e_iter)
//							>= /* *cfg::get().ds.RL * */2 && add) { //TODO: change to something reasonable.
////						missing_paired_info_count += TreatPairPairInfo<
////								graph_pack>(origin_gp, clustered_index,
////								first_info, second_info, 1);
//					}
//					if (origin_gp.g.length(*e_iter)
//							>= cfg::get().rr.max_repeat_length && clean) {
////						extra_paired_info_count +=
////								TreatPairPairInfo<graph_pack>(origin_gp,
////										clustered_index, first_info,
////										second_info, 0);
//					}
//				}
//			}
//		}
//	}
//	INFO(
//			"Paired info stats: missing = " << missing_paired_info_count
//					<< "; contradictional = " << extra_paired_info_count);
//}

}
