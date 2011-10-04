#ifndef IOPROCEDURES_HPP_
#define IOPROCEDURES_HPP_
#include <cmath>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <fstream>

#include "standard.hpp"
#include "logging.hpp"
#include "omni/paired_info.hpp"
#include "omni/omni_utils.hpp"

#include "omni/omni_tools.hpp"
#include "omni/omnigraph.hpp"

#include "omni/ID_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "omni/EdgeVertexFilter.hpp"
using namespace omnigraph;
using namespace debruijn_graph;

namespace omnigraph {
//todo think of inner namespace
//DECL_LOGGER("DataPrinter")

template<class Graph>
class DataPrinter {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	//	DataPrinter(/*const string& file_name,*/ Graph &g, IdTrackHandler<Graph> &old_IDs);
	void saveGraph(const string& file_name);
	void saveEdgeSequences(const string& file_name);
	void saveCoverage(const string& file_name);
	void savePaired(const string& file_name, PairedInfoIndex<Graph> const& PIIndex);
	void savePositions(const string& file_name,
			EdgesPositionHandler<Graph> const& EPHandler);

	void saveKmerMapper(const string& file_name,
			KmerMapper<K + 1, Graph> const& mapper);

	void close();

private:
	void save(FILE* file, EdgeId eid);
	//	void save(Sequence *sequence);
	void save(FILE* file, VertexId vid);
	Graph const& graph_;
	int edge_count_;
	//	map<EdgeId, typename IdTrackHandler<Graph>::realIdType> real_edge_ids_;
	IdTrackHandler<Graph> const& IdHandler_;
	EdgeVertexFilter<Graph> *filter_;
public:
	DataPrinter(/*const string& file_name,*/Graph const&g,
			IdTrackHandler<Graph> const& old_IDs) :
		graph_(g), IdHandler_(old_IDs) {
		DEBUG("Creating of saver started");
		edge_count_ = 0;
		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			edge_count_++;
		}
		filter_ = NULL;
	}
	DataPrinter(/*const string& file_name,*/Graph const& g,
			IdTrackHandler<Graph> const& old_IDs, EdgeVertexFilter<Graph> *filter) :
		graph_(g), IdHandler_(old_IDs), filter_(filter) {
		DEBUG("Creating of saver started");
		edge_count_ = 0;
		for (auto iter = filter_->EdgesBegin(); iter != filter_->EdgesEnd(); ++iter) {
			edge_count_++;
		}
	}
};

template<class Graph>
void DataPrinter<Graph>::saveGraph(const string& file_name) {

	FILE* file = fopen((file_name + ".grp").c_str(), "w");
	DEBUG("Graph saving to " << file_name << " started");
	VERIFY(file != NULL);
	if (filter_ == NULL) {
		int vertex_count = graph_.size();
		fprintf(file, "%d %d \n", vertex_count, edge_count_);
		for (auto iter = graph_.begin(); iter != graph_.end(); ++iter) {
			save(file, *iter);
		}

		fprintf(file, "\n");

		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			save(file, *iter);
		}
		DEBUG("Graph saving to " << file_name << " finished");
	} else {
		int vertex_count = filter_->VertexCount();
		fprintf(file, "%d %d \n", vertex_count, edge_count_);
		for (auto iter = filter_->VerticesBegin(); iter
				!= filter_->VerticesEnd(); ++iter) {
			save(file, *iter);
		}

		fprintf(file, "\n");

		for (auto iter = filter_->EdgesBegin(); iter != filter_->EdgesEnd(); ++iter) {
			save(file, *iter);
		}
		DEBUG("Graph saving to " << file_name << " finished");

	}
	fclose(file);
}

template<class Graph>
void DataPrinter<Graph>::save(FILE* file, VertexId vid) {
	fprintf(file, "%s\n", graph_.toPrint(vid, IdHandler_).c_str());
}

template<class Graph>
void DataPrinter<Graph>::save(FILE* file, EdgeId eid) {
	fprintf(file, "%s\n", graph_.toPrint(eid, IdHandler_).c_str());
}

template<class Graph>
void DataPrinter<Graph>::saveEdgeSequences(const string& file_name) {
	FILE* file = fopen((file_name + ".sqn").c_str(), "w");
	DEBUG("Saving sequences " << file_name <<" created");
	VERIFY(file != NULL);
	fprintf(file, "%d\n", edge_count_);
	if (filter_ == NULL) {
		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			fprintf(file, "%d ", IdHandler_.ReturnIntId(*iter));
			int len = graph_.EdgeNucls(*iter).size();
			for (int i = 0; i < len; i++)
				fprintf(file, "%c", nucl(graph_.EdgeNucls(*iter)[i]));
			fprintf(file, " .\n");
			//		fprintf(file, "%s .\n", graph_.EdgeNucls(*iter).str().c_str());
		}
	} else {
		for (auto iter = filter_->EdgesBegin(); iter != filter_->EdgesEnd(); ++iter) {
			fprintf(file, "%d ", IdHandler_.ReturnIntId(*iter));
			int len = graph_.EdgeNucls(*iter).size();
			for (int i = 0; i < len; i++)
				fprintf(file, "%c", nucl(graph_.EdgeNucls(*iter)[i]));
			fprintf(file, " .\n");
			//		fprintf(file, "%s .\n", graph_.EdgeNucls(*iter).str().c_str());
		}
	}
	fclose(file);
}

template<class Graph>
void DataPrinter<Graph>::saveCoverage(const string& file_name) {
	FILE* file = fopen((file_name + ".cvr").c_str(), "w");
	DEBUG("Saving coverage, " << file_name <<" created");
	VERIFY(file != NULL);
	fprintf(file, "%d\n", edge_count_);
	if (filter_ == NULL) {
		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			fprintf(file, "%d ", IdHandler_.ReturnIntId(*iter));
			fprintf(file, "%f .\n", graph_.coverage(*iter));
		}
	} else {
		for (auto iter = filter_->EdgesBegin(); iter != filter_->EdgesEnd(); ++iter) {
			fprintf(file, "%d ", IdHandler_.ReturnIntId(*iter));
			fprintf(file, "%f .\n", graph_.coverage(*iter));
		}
	}
	fclose(file);
}
/*
 template<class Graph>
 void DataPrinter<Graph>::saveIndex(const string& file_name) {
 FILE* file = fopen((file_name + ".ind").c_str(), "w");
 DEBUG("Saving index, " << file_name <<" created");
 VERIFY(file != NULL);
 fprintf(file, "%d\n", edge_count_);
 if (filter_ == NULL) {
 for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
 fprintf(file, "%d ", IdHandler_.ReturnIntId(*iter));
 fprintf(file, "%f .\n", graph_.coverage(*iter));
 }
 } else {
 for (auto iter = filter_->EdgesBegin(); iter != filter_->EdgesEnd(); ++iter) {
 fprintf(file, "%d ", IdHandler_.ReturnIntId(*iter));
 fprintf(file, "%f .\n", graph_.coverage(*iter));
 }
 }
 fclose(file);
 }
 */
template<class Graph>
void DataPrinter<Graph>::savePaired(const string& file_name,
		PairedInfoIndex<Graph> const& PIIndex) {
	FILE* file = fopen((file_name + ".prd").c_str(), "w");
	DEBUG("Saving paired info, " << file_name <<" created");
	VERIFY(file != NULL);
	if (filter_ == NULL) {
		fprintf(file, "%d\n", (int) PIIndex.size());
	} else {
		int filteredPIIsize = 0;
		for (auto iter = PIIndex.begin(); iter != PIIndex.end(); ++iter) {
			vector<PairInfo<typename Graph::EdgeId> > pair_infos = *iter;
			for (size_t i = 0; i < pair_infos.size(); i++) {
				if (filter_->EdgeIsPresent(pair_infos[i].first)
						&& filter_->EdgeIsPresent(pair_infos[i].second)) {
					filteredPIIsize++;
				}
			}
		}
		fprintf(file, "%d\n", filteredPIIsize);
	}
	for (auto iter = PIIndex.begin(); iter != PIIndex.end(); ++iter) {
		vector<PairInfo<typename Graph::EdgeId> > pair_infos = *iter;
		for (size_t i = 0; i < pair_infos.size(); i++) {
			if (filter_ == NULL) {
				fprintf(file, "%d %d %.2f %.2f %.2f .\n",
						IdHandler_.ReturnIntId(pair_infos[i].first),
						IdHandler_.ReturnIntId(pair_infos[i].second),
						pair_infos[i].d, pair_infos[i].weight,
						pair_infos[i].variance);
			} else {
				if (filter_->EdgeIsPresent(pair_infos[i].first)
						&& filter_->EdgeIsPresent(pair_infos[i].second)) {
					fprintf(file, "%d %d %.2f %.2f %.2f .\n",
							IdHandler_.ReturnIntId(pair_infos[i].first),
							IdHandler_.ReturnIntId(pair_infos[i].second),
							pair_infos[i].d, pair_infos[i].weight,
							pair_infos[i].variance);
				}
			}
		}
	}
	fclose(file);
}

template<class Graph>
void DataPrinter<Graph>::savePositions(const string& file_name,
		EdgesPositionHandler<Graph>const& EPHandler) {

	ofstream file((file_name + ".pos").c_str());

	DEBUG("Saving edges positions, " << file_name << " created");
	VERIFY(file != NULL);

	file << edge_count_ << endl;

	if (filter_ == NULL) {
		for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		    auto it = EPHandler.edges_positions().find(*iter);
		    VERIFY(it != EPHandler.edges_positions().end());

		    size_t size = it->second.size();
		    file << IdHandler_.ReturnIntId(*iter) << " " << size << endl;

			for (size_t i = 0; i < it->second.size(); i++)
				file << "    " <<it->second[i].contigId_<<": "<< it->second[i].start_ << " - " << it->second[i].end_ << endl;
		}
	} else {
		for (auto iter = filter_->EdgesBegin(); iter != filter_->EdgesEnd(); ++iter) {

		    auto it = EPHandler.edges_positions().find(*iter);
		    VERIFY(it != EPHandler.edges_positions().end());

			file << IdHandler_.ReturnIntId(*iter) << " " << it->second.size() << endl;

			for (size_t i = 0; i < it->second.size(); i++)
				file << "    "<<it->second[i].contigId_<<": " << it->second[i].start_ << " - " << it->second[i].end_ << endl;

		}
	}
}

template<class Graph>
void DataPrinter<Graph>::saveKmerMapper(const string& file_name,
		KmerMapper<K + 1, Graph> const& mapper) {

	std::ofstream file;
	file.open((file_name + ".kmm").c_str(),  std::ios_base::binary | std::ios_base::out);
	DEBUG("Saving kmer mapper, " << file_name <<" created");
	VERIFY(file.is_open());

	u_int32_t k_ = K;
	file.write((char *) &k_, sizeof(u_int32_t));
	mapper.BinWrite(file);

	file.close();
}

template<class Graph>
class DataScanner {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	//	DataPrinter(/*const string& file_name,*/ Graph &g, IdTrackHandler<Graph> &old_IDs);
	void loadNonConjugateGraph(const string& file_name, bool with_Sequence);
	void loadConjugateGraph(const string& file_name, bool with_Sequence);

	//	void saveEdgeSequences(const string& file_name);
	void loadCoverage(const string& file_name);
	void loadPaired(const string& file_name, PairedInfoIndex<Graph>& PIIndex);
	void loadPositions(const string& file_name,
			EdgesPositionHandler<Graph>& EPHandler);

	void loadKmerMapper(const string& file_name,
			KmerMapper<K + 1, Graph>& mapper);

	void close();

private:
	//

	Graph &graph_;
	int edge_count_;
	//	map<EdgeId, typename IdTrackHandler<Graph>::realIdType> real_edge_ids_;
	IdTrackHandler<Graph>& IdHandler_;
public:
	DataScanner(/*const string& file_name,*/Graph &g, IdTrackHandler<Graph>&new_IDs)
	    : graph_    (g)
	    , IdHandler_(new_IDs)
	{
		INFO("Creating of scanner started");
		edge_count_ = 0;
	}
};

template<class Graph>
void DataScanner<Graph>::loadNonConjugateGraph(const string& file_name,
		bool with_Sequence) {
	int read_count;
	FILE* file = fopen((file_name + ".grp").c_str(), "r");
	if (file == NULL) WARN("File "<<(file_name + ".grp")<<" not found");
	VERIFY(file != NULL);
	FILE* sequence_file = fopen((file_name + ".sqn").c_str(), "r");
	VERIFY(sequence_file != NULL);

	INFO("Reading NON conjugate de bruujn graph from " << file_name << " started");
	int vertex_count;
	read_count = fscanf(file, "%d %d \n", &vertex_count, &edge_count_);
	VERIFY(read_count == 2);
	for (int i = 0; i < vertex_count; i++) {
		int vertex_real_id;
		read_count = fscanf(file, "Vertex %d", &vertex_real_id);
		VERIFY(read_count == 1);
		char c = 'a';
		while (c != '.') {
			read_count = fscanf(file, "%c", &c);
			VERIFY(read_count == 1);
		}
		read_count = fscanf(file, "\n");
		VERIFY(read_count == 0);
		VertexId vid = graph_.AddVertex();
		IdHandler_.AddVertexIntId(vid, vertex_real_id);
		TRACE(vid);
	}
	int tmp_edge_count;
	read_count = fscanf(sequence_file, "%d", &tmp_edge_count);
	VERIFY(read_count == 1);
	VERIFY(edge_count_ == tmp_edge_count);
	char longstring[1000500];
	for (int i = 0; i < edge_count_; i++) {
		int e_real_id, start_id, fin_id, length;
		read_count = fscanf(file, "Edge %d : %d -> %d, l = %d", &e_real_id,
				&start_id, &fin_id, &length);
		VERIFY(read_count == 4);
		read_count = fscanf(sequence_file, "%d %s .", &e_real_id, longstring);
		VERIFY(read_count == 2);
		//does'nt matter, whether it was conjugate or not.
		char c = 'a';
		while (c != '.') {
			read_count = fscanf(file, "%c", &c);
			VERIFY(read_count == 1);
		}
		read_count = fscanf(file, "\n");
		VERIFY(read_count == 0);
		Sequence tmp(longstring);
		TRACE(start_id<<" "<< fin_id <<" "<< IdHandler_.ReturnVertexId(start_id)<<" "<< IdHandler_.ReturnVertexId(fin_id));
		EdgeId eid = graph_.AddEdge(IdHandler_.ReturnVertexId(start_id),
				IdHandler_.ReturnVertexId(fin_id), tmp);
		IdHandler_.AddEdgeIntId(eid, e_real_id);

	}
	fclose(file);
	fclose(sequence_file);
}

template<class Graph>
void DataScanner<Graph>::loadConjugateGraph(const string& file_name,
		bool with_Sequence) {
	int read_count;
	FILE* file = fopen((file_name + ".grp").c_str(), "r");
	VERIFY(file != NULL);
	FILE* sequence_file = fopen((file_name + ".sqn").c_str(), "r");
	VERIFY(sequence_file != NULL);
	set<int> vertex_set;
	set<int> edge_set;
	INFO("Reading conjugate de bruijn  graph from " << file_name << " started");
	int vertex_count;
	read_count = fscanf(file, "%d %d \n", &vertex_count, &edge_count_);
	VERIFY(read_count == 2);
	for (int i = 0; i < vertex_count; i++) {
		int vertex_real_id, conjugate_id;
		read_count = fscanf(file, "Vertex %d ~ %d .\n", &vertex_real_id,
				&conjugate_id);
		TRACE("Vertex "<<vertex_real_id<<" ~ "<<conjugate_id<<" .");
		VERIFY(read_count == 2);

		if (vertex_set.find(vertex_real_id) == vertex_set.end()) {
			VertexId vid = graph_.AddVertex();
			VertexId conj_vid = graph_.conjugate(vid);

			IdHandler_.AddVertexIntId(vid, vertex_real_id);
			IdHandler_.AddVertexIntId(conj_vid, conjugate_id);
			vertex_set.insert(conjugate_id);
			TRACE(vid<<" ( "<< IdHandler_.ReturnVertexId(vertex_real_id) <<" )   "<< conj_vid << "( "<<IdHandler_.ReturnVertexId(conjugate_id)<<" )  added");
		}
	}
	int tmp_edge_count;
	read_count = fscanf(sequence_file, "%d", &tmp_edge_count);
	VERIFY(read_count == 1);
	VERIFY(edge_count_ == tmp_edge_count);
	char longstring[1000500];
	for (int i = 0; i < edge_count_; i++) {
		int e_real_id, start_id, fin_id, length, conjugate_edge_id;
		read_count = fscanf(file, "Edge %d : %d -> %d, l = %d ~ %d .\n",
				&e_real_id, &start_id, &fin_id, &length, &conjugate_edge_id);
		VERIFY(read_count == 5);
		read_count = fscanf(sequence_file, "%d %s .", &e_real_id, longstring);
		VERIFY(read_count == 2);
		TRACE("Edge "<<e_real_id<<" : "<<start_id<<" -> " << fin_id << " l = " << length << " ~ "<< conjugate_edge_id);
		if (edge_set.find(e_real_id) == edge_set.end()) {
			Sequence tmp(longstring);
			TRACE(start_id<<" "<< fin_id <<" "<< IdHandler_.ReturnVertexId(start_id)<<" "<< IdHandler_.ReturnVertexId(fin_id));
			EdgeId eid = graph_.AddEdge(IdHandler_.ReturnVertexId(start_id),
					IdHandler_.ReturnVertexId(fin_id), tmp);
			IdHandler_.AddEdgeIntId(eid, e_real_id);
			IdHandler_.AddEdgeIntId(graph_.conjugate(eid), conjugate_edge_id);
			edge_set.insert(conjugate_edge_id);

		}

	}
	fclose(file);
	fclose(sequence_file);
}

template<class Graph>
void DataScanner<Graph>::loadCoverage(const string& file_name) {
	int read_count;
	FILE* file = fopen((file_name + ".cvr").c_str(), "r");
	VERIFY(file != NULL);
	INFO("Reading coverage from " << file_name << " started");
	int edge_count;
	read_count = fscanf(file, "%d \n", &edge_count);
	VERIFY(read_count == 1);
	VERIFY(edge_count == edge_count_);
	for (int i = 0; i < edge_count; i++) {
		int edge_real_id;
		double edge_coverage;
		read_count = fscanf(file, "%d %lf .\n", &edge_real_id, &edge_coverage);
		VERIFY(read_count == 2);
		TRACE(edge_real_id<< " "<<edge_coverage <<" . ");
		EdgeId eid = IdHandler_.ReturnEdgeId(edge_real_id);
		TRACE("EdgeId "<<eid);
		graph_.coverage_index().SetCoverage(eid, edge_coverage * graph_.length(eid));
	}
	fclose(file);
}

template<class Graph>
void DataScanner<Graph>::loadPaired(const string& file_name,
		PairedInfoIndex<Graph>& PIIndex) {
	int read_count;
	FILE* file = fopen((file_name + ".prd").c_str(), "r");
	DEBUG((file_name + ".prd"));
	VERIFY(file != NULL);
	INFO("Reading paired info from " << file_name << " started");
	int paired_count;
	read_count = fscanf(file, "%d \n", &paired_count);
	VERIFY(read_count == 1);
	for (int i = 0; i < paired_count; i++) {
		int first_real_id, second_real_id;
		double w, d, v;
		read_count = fscanf(file, "%d %d %lf %lf %lf .\n", &first_real_id,
				&second_real_id, &d, &w, &v);
		VERIFY(read_count == 5);
		TRACE(first_real_id<< " " << second_real_id << " " << d << " " << w << " " << v);
		TRACE (IdHandler_.ReturnEdgeId(first_real_id)<<" "<< IdHandler_.ReturnEdgeId(second_real_id)<<" "<< d<<" "<< w);
		PairInfo<typename Graph::EdgeId> *p_info = new PairInfo<
				typename Graph::EdgeId> (
				IdHandler_.ReturnEdgeId(first_real_id),
				IdHandler_.ReturnEdgeId(second_real_id), d, w, v);
		PIIndex.AddPairInfo(*p_info, 0);
	}
	DEBUG("PII SIZE " << PIIndex.size());
	fclose(file);
}

template<class Graph>
void DataScanner<Graph>::loadPositions(const string& file_name,
		EdgesPositionHandler<Graph>& EPHandler) {
	int read_count;
	FILE* file = fopen((file_name + ".pos").c_str(), "r");
	VERIFY(file != NULL);
	DEBUG("Reading edges positions, " << file_name <<" started");
	VERIFY(file != NULL);
	int pos_count;
	read_count = fscanf(file, "%d\n", &pos_count);
	VERIFY(read_count == 1);
	for (int i = 0; i < pos_count; i++) {
		int edge_real_id, pos_info_count, contigId;
		read_count = fscanf(file, "%d %d\n", &edge_real_id, &pos_info_count);
		VERIFY(read_count == 2);
		for (int j = 0; j < pos_info_count; j++) {
			int start_pos, end_pos;
			read_count = fscanf(file, "%d: %d - %d \n", &contigId, &start_pos, &end_pos);
			VERIFY(read_count == 3);
			EdgeId eid = IdHandler_.ReturnEdgeId(edge_real_id);
			EPHandler.AddEdgePosition(eid, start_pos, end_pos, contigId);
		}
	}
	fclose(file);
}

template<class Graph>
void DataScanner<Graph>::loadKmerMapper(const string& file_name,
		KmerMapper<K + 1, Graph>& mapper) {

	mapper.clear();
	std::ifstream file;
	file.open((file_name + ".kmm").c_str(), std::ios_base::binary | std::ios_base::in);
	DEBUG("Reading kmer mapper, " << file_name <<" started");
	VERIFY(file.is_open());

	u_int32_t k_;
	file.read((char *) &k_, sizeof(u_int32_t));

	VERIFY_MSG(k_ == K, "Cannot read kmer mapper, different Ks");
	mapper.BinRead(file);

	file.close();
}




template<class Graph>
void printGraph(Graph const& g, IdTrackHandler<Graph> &old_IDs,
		const string &file_name, PairedInfoIndex<Graph> &paired_index,
		EdgesPositionHandler<Graph> &edges_positions,
		EdgeVertexFilter<Graph> *filter) {
	DataPrinter<Graph> dataPrinter(g, old_IDs, filter);
	dataPrinter.saveGraph(file_name);
	dataPrinter.saveEdgeSequences(file_name);
	dataPrinter.saveCoverage(file_name);
	dataPrinter.savePaired(file_name, paired_index);
	dataPrinter.savePositions(file_name, edges_positions);
}

template<class Graph>
void printGraph(Graph const & g, IdTrackHandler<Graph> const& old_IDs,
		const string &file_name, PairedInfoIndex<Graph> const& paired_index,
		EdgesPositionHandler<Graph> const& edges_positions,
		PairedInfoIndex<Graph> const* etalon_index = 0,
		PairedInfoIndex<Graph> const* clustered_index = 0,
		KmerMapper<K + 1, Graph> const* mapper = 0) {

	DataPrinter<Graph> dataPrinter(g, old_IDs);
	dataPrinter.saveGraph(file_name);
	dataPrinter.saveEdgeSequences(file_name);
	dataPrinter.saveCoverage(file_name);
	dataPrinter.savePaired(file_name, paired_index);
	//todo delete
	if (etalon_index) {
		dataPrinter.savePaired(file_name + "_et", *etalon_index);
	}
	if (clustered_index) {
		dataPrinter.savePaired(file_name + "_cl", *clustered_index);
	}
	dataPrinter.savePositions(file_name, edges_positions);

	if (mapper) {
		dataPrinter.saveKmerMapper(file_name, *mapper);
	}
}

template<class Graph>
void printGraph(Graph const & g, IdTrackHandler<Graph> &old_IDs,
		const string &file_name, PairedInfoIndex<Graph> &paired_index) {
	DataPrinter<Graph> dataPrinter(g, old_IDs);
	dataPrinter.saveGraph(file_name);
	dataPrinter.saveEdgeSequences(file_name);
	dataPrinter.saveCoverage(file_name);
	dataPrinter.savePaired(file_name, paired_index);

}

template<class Graph>
void printKmerMapper(Graph const & g, IdTrackHandler<Graph> const& old_IDs,
		const string &file_name, KmerMapper<K + 1, Graph>& mapper) {

	DataPrinter<Graph> dataPrinter(g, old_IDs);
	dataPrinter.saveKmerMapper(file_name, mapper);
}


template<class Graph>
void scanNCGraph(Graph & g, IdTrackHandler<Graph>&new_IDs,
		const string &file_name, PairedInfoIndex<Graph>* paired_index,
		EdgesPositionHandler<Graph> &edges_positions,
		PairedInfoIndex<Graph>* etalon_index = 0,
		PairedInfoIndex<Graph>* clustered_index = 0) {
	DataScanner<Graph> dataScanner(g, new_IDs);
	dataScanner.loadNonConjugateGraph(file_name, true);
	dataScanner.loadCoverage(file_name);
	if (paired_index) {
		dataScanner.loadPaired(file_name, *paired_index);
	}
	if (etalon_index) {
		dataScanner.loadPaired(file_name + "_et", *etalon_index);
	}
	if (clustered_index) {
		dataScanner.loadPaired(file_name + "_cl", *clustered_index);
	}
	dataScanner.loadPositions(file_name, edges_positions);
}

template<class Graph>
void scanNCGraph(Graph & g, IdTrackHandler<Graph> &new_IDs,
		const string &file_name, PairedInfoIndex<Graph>& paired_index) {
	DataScanner<Graph> dataScanner(g, new_IDs);
	dataScanner.loadNonConjugateGraph(file_name, true);
	dataScanner.loadCoverage(file_name);
	dataScanner.loadPaired(file_name, paired_index);
}

template<class Graph>
void scanConjugateGraph(Graph * g, IdTrackHandler<Graph> *new_IDs,
		const string &file_name, PairedInfoIndex<Graph>* paired_index = 0,
		EdgesPositionHandler<Graph> *edges_positions = NULL,
		PairedInfoIndex<Graph>* etalon_index = 0,
		PairedInfoIndex<Graph>* clustered_index = 0,
		KmerMapper<K + 1, Graph> * mapper = 0) {
	//ToDo Apply * vs & conventions
	DataScanner<Graph> dataScanner(*g, *new_IDs);
	dataScanner.loadConjugateGraph(file_name, true);
	dataScanner.loadCoverage(file_name);
	if (paired_index) {
		dataScanner.loadPaired(file_name, *paired_index);
	}
	if (edges_positions != NULL)
		dataScanner.loadPositions(file_name, *edges_positions);
	if (etalon_index) {
		dataScanner.loadPaired(file_name + "_et", *etalon_index);
	}
	if (clustered_index) {
		dataScanner.loadPaired(file_name + "_cl", *clustered_index);
	}
	if (mapper) {
		dataScanner.loadKmerMapper(file_name, *mapper);
	}
}

template<class Graph>
void scanKmerMapper(Graph& g, IdTrackHandler<Graph>& new_IDs,
		const string &file_name, KmerMapper<K + 1, Graph> * mapper) {

	DataScanner<Graph> dataScanner(g, new_IDs);
	dataScanner.loadKmerMapper(file_name, *mapper);
}


}
#endif /* IOPROCEDURES_HPP_ */
