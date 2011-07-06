#ifndef IOPROCEDURES_HPP_
#define IOPROCEDURES_HPP_
#include <cmath>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "logging.hpp"
#include "paired_info.hpp"
#include "omni_utils.hpp"

#include "omni_tools.hpp"
#include "omnigraph.hpp"

#include "debruijn/ID_track_handler.hpp"
using namespace omnigraph;
using namespace debruijn_graph;

namespace omnigraph {
//DECL_LOGGER("DataPrinter")

template <class Graph>
class DataPrinter {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
//	DataPrinter(/*const string& file_name,*/ Graph &g, IdTrackHandler<Graph> &old_IDs);
	void saveGraph(const string& file_name);
	void saveEdgeSequences(const string& file_name);
	void saveCoverage(const string& file_name);
	void savePaired(const string& file_name, PairedInfoIndex<Graph>& PIIndex);
	void close();

private:
	void save(FILE* file, EdgeId eid);
//	void save(Sequence *sequence);
	void save(FILE* file, VertexId vid);

	Graph &graph_;
	int edge_count_;
//	map<EdgeId, typename IdTrackHandler<Graph>::realIdType> real_edge_ids_;
	IdTrackHandler<Graph> &IdHandler_;
public:
	DataPrinter(/*const string& file_name,*/ Graph &g, IdTrackHandler<Graph> &old_IDs) : graph_(g), IdHandler_(old_IDs) {
		INFO("Creating of saver started");
		edge_count_ = 0;
		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			edge_count_++;
		}
	}
};

template<class Graph>
void DataPrinter<Graph>::saveGraph(const string& file_name) {

	FILE* file = fopen((file_name + ".grp").c_str(), "w");
	INFO("Graph saving to " << file_name << " started");
	assert(file != NULL);
	int vertex_count = graph_.size();
	fprintf(file, "%d %d \n", vertex_count, edge_count_);
	for (auto iter = graph_.begin(); iter != graph_.end(); ++iter) {
		save(file, *iter);
	}

	fprintf(file, "\n");

	for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		save(file, *iter);
	}
	INFO("Graph saving to " << file_name << " finished");

	fclose(file);
}

template<class Graph>
void DataPrinter<Graph>::save(FILE* file, VertexId vid) {
	fprintf(file, "%s\n", graph_.toPrint(vid, IdHandler_).c_str());
}

template<class Graph>
void DataPrinter<Graph>::save(FILE* file, EdgeId eid) {
	fprintf(file, "%s\n",  graph_.toPrint(eid, IdHandler_).c_str());
}

template<class Graph>
void DataPrinter<Graph>::saveEdgeSequences(const string& file_name) {
	FILE* file = fopen((file_name + ".sqn").c_str(), "w");
	DEBUG("Saving sequences " << file_name <<" created");
	assert(file != NULL);
	fprintf(file, "%d\n", edge_count_);
	for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		fprintf(file, "%d ", IdHandler_.ReturnIntId(*iter));
		int len = graph_.EdgeNucls(*iter).size();
		for(int i = 0; i < len; i++ )
			fprintf(file, "%c",nucl(graph_.EdgeNucls(*iter)[i]));
		fprintf(file, " .\n");
//		fprintf(file, "%s .\n", graph_.EdgeNucls(*iter).str().c_str());
	}
	fclose(file);
}

template<class Graph>
void DataPrinter<Graph>::saveCoverage(const string& file_name) {
	FILE* file = fopen((file_name+".cvr").c_str(), "w");
	DEBUG("Saving coverage, " << file_name <<" created");
	assert(file != NULL);
	fprintf(file, "%d\n", edge_count_);
	for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		fprintf(file, "%d ", IdHandler_.ReturnIntId(*iter));
		fprintf(file, "%f .\n", graph_.coverage(*iter));
	}
	fclose(file);
}

template<class Graph>
void DataPrinter<Graph>::savePaired(const string& file_name, PairedInfoIndex<Graph>& PIIndex) {
	FILE* file = fopen((file_name + ".prd").c_str(), "w");
	DEBUG("Saving paired info, " << file_name <<" created");
	assert(file != NULL);
	fprintf(file, "%d\n", PIIndex.size());
	for (auto iter = PIIndex.begin(); iter != PIIndex.end(); ++iter) {
		vector<PairInfo<typename Graph::EdgeId> > pair_infos = *iter;
		for(size_t i = 0; i < pair_infos.size(); i++)
			fprintf(file, "%d %d %.0f %.0f .\n", IdHandler_.ReturnIntId(pair_infos[i].first), IdHandler_.ReturnIntId(pair_infos[i].second), pair_infos[i].d, pair_infos[i].weight);
	}
	fclose(file);
}

template <class Graph>
class DataScanner {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
//	DataPrinter(/*const string& file_name,*/ Graph &g, IdTrackHandler<Graph> &old_IDs);
	void loadNonConjugateGraph(const string& file_name, bool with_Sequence);
//	void saveEdgeSequences(const string& file_name);
	void loadCoverage(const string& file_name);
	void loadPaired(const string& file_name, PairedInfoIndex<Graph>& PIIndex);
	void close();

private:
//

	Graph &graph_;
	int edge_count_;
//	map<EdgeId, typename IdTrackHandler<Graph>::realIdType> real_edge_ids_;
	IdTrackHandler<Graph> &IdHandler_;
public:
	DataScanner(/*const string& file_name,*/ Graph &g, IdTrackHandler<Graph> &new_IDs) : graph_(g), IdHandler_(new_IDs) {
		INFO("Creating of scanner started");
		edge_count_ = 0;
	}
};
template<class Graph>
void DataScanner<Graph>::loadNonConjugateGraph(const string& file_name, bool with_Sequence) {

	FILE* file = fopen((file_name + ".grp").c_str(), "r");
	assert(file != NULL);
	FILE* sequence_file = fopen((file_name + ".sqn").c_str(), "r");
	assert(sequence_file != NULL);

	INFO("Reading graph from " << file_name << " started");
	int vertex_count;
	assert(fscanf(file, "%d %d \n", &vertex_count, &edge_count_) == 2);
	for (int i = 0; i < vertex_count; i++) {
		int vertex_real_id;
		assert(fscanf(file, "Vertex %d", &vertex_real_id) == 1);
		char c = 'a';
		while (c != '.')
			assert(fscanf(file, "%c", &c) == 1);
		assert( fscanf(file, "\n") == 0);
		VertexId vid = graph_.AddVertex();
		IdHandler_.AddVertexIntId(vid, vertex_real_id);
		DEBUG(vid);
	}
	int tmp_edge_count;
	assert(fscanf(sequence_file, "%d", &tmp_edge_count) == 1);
	assert(edge_count_ == tmp_edge_count);
	char longstring[1000500];
	for (int i = 0; i < edge_count_; i++){
		int e_real_id, start_id, fin_id, length;
		assert(fscanf(file, "Edge %d : %d -> %d, l = %d", &e_real_id, &start_id, &fin_id, &length) == 4);
		assert(fscanf(sequence_file, "%d %s .", &e_real_id, longstring) == 2);

		//does'nt matter, whether it was conjugate or not.
		char c = 'a';
		while (c != '.')
			assert(fscanf(file, "%c", &c) == 1);
		assert( fscanf(file, "\n") == 0);
		Sequence tmp(longstring);
		DEBUG(start_id<<" "<<  fin_id <<" "<< IdHandler_.ReturnVertexId(start_id)<<" "<< IdHandler_.ReturnVertexId(fin_id));
		EdgeId eid = graph_.AddEdge(IdHandler_.ReturnVertexId(start_id), IdHandler_.ReturnVertexId(fin_id), tmp);
		IdHandler_.AddEdgeIntId(eid, e_real_id);

	}


	fclose(file);
	fclose(sequence_file);

}

template<class Graph>
void DataScanner<Graph>::loadCoverage(const string& file_name) {

	FILE* file = fopen((file_name + ".cvr").c_str(), "r");
	assert(file != NULL);
	INFO("Reading coverage from " << file_name << " started");
	int edge_count;
	assert(fscanf(file, "%d \n", &edge_count) == 1);
	assert(edge_count == edge_count_);
	for (int i = 0; i < edge_count; i++) {
		int edge_real_id;
		double edge_coverage;
		assert(fscanf(file, "%d %lf .\n", &edge_real_id, &edge_coverage) == 2);
		EdgeId eid = IdHandler_.ReturnEdgeId(edge_real_id);
		graph_.SetCoverage(eid, edge_coverage * graph_.length(eid));
	}
	fclose(file);
}

template<class Graph>
void DataScanner<Graph>::loadPaired(const string& file_name, PairedInfoIndex<Graph>& PIIndex) {

	FILE* file = fopen((file_name + ".prd").c_str(), "r");
	assert(file != NULL);
	INFO("Reading paired info from " << file_name << " started");
	int paired_count;
	assert(fscanf(file, "%d \n", &paired_count) == 1);
	for (int i = 0; i < paired_count; i++) {
		int first_real_id, second_real_id;
		double w, d;
		assert(fscanf(file, "%d %d %lf %lf .\n", &first_real_id, &second_real_id, &d, &w) == 4);
		DEBUG(first_real_id<< " " << second_real_id << " " << d << " " << w);
		DEBUG (IdHandler_.ReturnEdgeId(first_real_id)<<" "<< IdHandler_.ReturnEdgeId(second_real_id)<<" "<< d<<" "<< w);
		PairInfo<typename Graph::EdgeId> *p_info = new PairInfo<typename Graph::EdgeId>(IdHandler_.ReturnEdgeId(first_real_id), IdHandler_.ReturnEdgeId(second_real_id), d, w);
		PIIndex.AddPairInfo(*p_info, 0);
	}
	DEBUG("PII SIZE " << PIIndex.size());
	fclose(file);
}

}
#endif /* IOPROCEDURES_HPP_ */
