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

	FILE* file = fopen(file_name.c_str(), "w");
	INFO("Graph saving to " << file_name << "started");
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
	INFO("Graph saving to " << file_name << "finished");

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
	FILE* file = fopen(file_name.c_str(), "w");
	DEBUG("Saving sequences " << file_name <<" created");
	assert(file != NULL);
	fprintf(file, "%d\n", edge_count_);
	for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		fprintf(file, "%d ", IdHandler_.ReturnIntId(*iter));
	}
	fclose(file);
}
}
/*
class DataReader {
	FILE *f_;
public:
	DataReader(char *fileName);
	DataReader(const char *fileName);
	int read(int &a);
	int read(long long &a);
	int read(Edge * &edge);
	int read(Sequence * &sequence);
	void read(VertexPrototype * &v);
	void readLongEdgesMap(longEdgesMap &map);
	void readIntArray(int *array, int length);
	void readIntArray(int *array, int length, int width);
	template<typename keyType, typename valueType>
	void read(map<keyType, valueType> &m);
	template<typename valueType>
	void read(vector<valueType> &v);
	void close();
};

template<typename keyType, typename valueType>
void DataPrinter::output(map<keyType, valueType> m) {
	output((int)m.size());
	for(typename map<keyType, valueType>::iterator it = m.begin(); it != m.end(); ++it) {
		output(it->first);
		output(it->second);
	}
	fprintf(f_, "\n");
}

template<typename keyType, typename valueType>
void DataReader::read(map<keyType, valueType> &m) {
	m.clear();
	int size;
	read(size);
	keyType key;
	valueType value;
	for(int i = 0; i < size; i++) {
		read(key);
		read(value);
		m[key] = value;
	}
	fscanf(f_, "\n");
}

template<typename valueType>
void DataPrinter::output(vector<valueType> v) {
	output((int)v.size());
	for(typename vector<valueType>::iterator it = v.begin(); it != v.end(); ++it) {
		output(*it);
	}
	fprintf(f_, "\n");
}

template<typename valueType>
void DataReader::read(vector<valueType> &v) {
	v.clear();
	int size;
	read(size);
	v.reserve(size);
	valueType value;
	for(int i = 0; i < size; i++) {
		read(value);
		v.push_back(value);
	}
	fscanf(f_, "\n");
}

template<class Graph>
void save(char *fileName, Graph &g);
template<class Graph>
void load(char *fileName, Graph &g);

template<class Graph>
void save(string fileName, Graph &g);
template<class Graph>
void load(string fileName, Graph &g);
*/
#endif /* IOPROCEDURES_HPP_ */
