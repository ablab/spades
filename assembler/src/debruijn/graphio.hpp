//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <cmath>
#include <set>
#include <map>
#include <algorithm>
#include <fstream>
#include <cstdio>

#include "standard.hpp"
#include "omni/paired_info.hpp"
#include "omni/omni_utils.hpp"
#include "omni/abstract_conjugate_graph.hpp"
#include "omni/abstract_nonconjugate_graph.hpp"
#include "utils.hpp"
#include "new_debruijn.hpp"

#include "omni/omni_tools.hpp"
#include "omni/omnigraph.hpp"

#include "omni/id_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "omni/graph_component.hpp"

namespace debruijn_graph {
using namespace omnigraph;
//todo think of inner namespace

template<class KmerMapper>
void SaveKmerMapper(const string& file_name,
		const KmerMapper& mapper) {
	std::ofstream file;
	file.open((file_name + ".kmm").c_str(),
			std::ios_base::binary | std::ios_base::out);
	DEBUG("Saving kmer mapper, " << file_name <<" created");
	VERIFY(file.is_open());

	u_int32_t k_ = KmerMapper::k_value;
	file.write((char *) &k_, sizeof(u_int32_t));
	mapper.BinWrite(file);

	file.close();
}

template<class KmerMapper>
void LoadKmerMapper(const string& file_name,
		KmerMapper& kmer_mapper) {
	kmer_mapper.clear();
	std::ifstream file;
	file.open((file_name + ".kmm").c_str(),
			std::ios_base::binary | std::ios_base::in);
	DEBUG("Reading kmer mapper, " << file_name <<" started");
	VERIFY(file.is_open());

	u_int32_t k_;
	file.read((char *) &k_, sizeof(u_int32_t));

	VERIFY_MSG(k_ == KmerMapper::k_value, "Cannot read kmer mapper, different Ks");
	kmer_mapper.BinRead(file);

	file.close();
}

template<class Graph>
class DataPrinter {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	void saveGraph(const string& file_name);
	void saveEdgeSequences(const string& file_name);
	void saveCoverage(const string& file_name);
	void savePaired(const string& file_name,
			PairedInfoIndex<Graph> const& paired_index);
	void savePositions(const string& file_name,
			EdgesPositionHandler<Graph> const& ref_pos);

private:
	void save(FILE* file, EdgeId eid);
	void save(FILE* file, VertexId vid);

	const GraphComponent<Graph> component_;
	const BaseIdTrackHandler<VertexId, EdgeId>& int_ids_;

	virtual std::string toPrint(VertexId v) const = 0;
	virtual std::string toPrint(EdgeId e) const = 0;

protected:
//	DataPrinter(Graph const& g, IdTrackHandler<Graph> const& int_ids) :
//			component_(g), int_ids_(int_ids) {
//		DEBUG("Creating of saver started");
//		edge_count_ = 0;
//		if (graph_component) {
//			for (auto iter = graph_component_->EdgesBegin();
//					iter != graph_component_->EdgesEnd(); ++iter) {
//				edge_count_++;
//			}
//		} else {
//			for (auto iter = graph_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
//				edge_count_++;
//			}
//		}
//	}

//todo optimize component copy
	DataPrinter(const GraphComponent<Graph>& component,
			BaseIdTrackHandler<VertexId, EdgeId> const& int_ids) :
			component_(component), int_ids_(int_ids) {
	}

//	template<class VertexIt>
//	DataPrinter(const Graph& g, VertexIt begin, VertexIt end,
//			IdTrackHandler<Graph> const& int_ids, bool conjugate) :
//			component_(g, begin, end, conjugate), int_ids_(int_ids) {
//	}

	const GraphComponent<Graph>& component() const {
		return component_;
	}

	const BaseIdTrackHandler<VertexId, EdgeId>& id_handler() const {
		return int_ids_;
	}

public:

	virtual ~DataPrinter() {

	}
};

template<class Graph>
void DataPrinter<Graph>::saveGraph(const string& file_name) {

	FILE* file = fopen((file_name + ".grp").c_str(), "w");
	DEBUG("Graph saving to " << file_name << " started");
	VERIFY_MSG(file != NULL,
	"Couldn't open file " << (file_name + ".grp") << " on write");
	size_t vertex_count = component_.v_size();
	size_t edge_count = component_.e_size();
	fprintf(file, "%ld %ld \n", vertex_count, edge_count);
	for (auto iter = component_.v_begin(); iter != component_.v_end(); ++iter) {
		save(file, *iter);
	}

	fprintf(file, "\n");

	for (auto iter = component_.e_begin(); iter != component_.e_end(); ++iter) {
		save(file, *iter);
	}DEBUG("Graph saving to " << file_name << " finished");

	fclose(file);
}

template<class Graph>
void DataPrinter<Graph>::save(FILE* file, VertexId vid) {
	fprintf(file, "%s\n", toPrint(vid).c_str());
}

template<class Graph>
void DataPrinter<Graph>::save(FILE* file, EdgeId eid) {
	fprintf(file, "%s\n", toPrint(eid).c_str());
}

template<class Graph>
void DataPrinter<Graph>::saveEdgeSequences(const string& file_name) {
	FILE* file = fopen((file_name + ".sqn").c_str(), "w");
	DEBUG("Saving sequences " << file_name <<" created");
	VERIFY(file != NULL);
	//fprintf(file, "%ld\n", component_.e_size());
	for (auto iter = component_.e_begin(); iter != component_.e_end(); ++iter) {
		fprintf(file, ">%d\n", int_ids_.ReturnIntId(*iter));
		int len = component_.g().EdgeNucls(*iter).size();
		for (int i = 0; i < len; i++)
			fprintf(file, "%c", nucl(component_.g().EdgeNucls(*iter)[i]));
		fprintf(file, "\n");
		//		fprintf(file, "%s .\n", graph_.EdgeNucls(*iter).str().c_str());
	}
	fclose(file);
}

template<class Graph>
void DataPrinter<Graph>::saveCoverage(const string& file_name) {
	FILE* file = fopen((file_name + ".cvr").c_str(), "w");
	DEBUG("Saving coverage, " << file_name <<" created");
	VERIFY(file != NULL);
	fprintf(file, "%ld\n", component_.e_size());
	for (auto iter = component_.e_begin(); iter != component_.e_end(); ++iter) {
		fprintf(file, "%d ", int_ids_.ReturnIntId(*iter));
		fprintf(file, "%f .\n", component_.g().coverage(*iter));
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
		PairedInfoIndex<Graph> const& paired_index) {
	typedef vector<PairInfo<typename Graph::EdgeId>> Infos;
	FILE* file = fopen((file_name + ".prd").c_str(), "w");
	DEBUG("Saving paired info, " << file_name <<" created");
	VERIFY(file != NULL);

	Infos to_save;
	for (auto it = component_.e_begin(); it != component_.e_end(); ++it) {
		Infos infos = paired_index.GetEdgeInfo(*it);
		for (auto info_it = infos.begin(); info_it != infos.end(); ++info_it) {
			if (info_it->d < 0) { continue; }
			if (component_.contains(info_it->second)) {
				to_save.push_back(*info_it);
			}
		}
	}

	fprintf(file, "%ld\n", to_save.size());

	for (auto it = to_save.begin(); it != to_save.end(); ++it) {
		fprintf(file, "%d %d %.2f %.2f %.2f .\n",
				int_ids_.ReturnIntId(it->first),
				int_ids_.ReturnIntId(it->second), it->d, it->weight,
				it->variance);
	}

	fclose(file);
}

template<class Graph>
void DataPrinter<Graph>::savePositions(const string& file_name,
		EdgesPositionHandler<Graph> const& ref_pos) {

	ofstream file((file_name + ".pos").c_str());

	DEBUG("Saving edges positions, " << file_name << " created");
	VERIFY(file != NULL);

	file << component_.e_size() << endl;

	for (auto it = component_.e_begin(); it != component_.e_end(); ++it) {

		auto pos_it = ref_pos.edges_positions().find(*it);
		VERIFY(pos_it != ref_pos.edges_positions().end());

		file << id_handler().ReturnIntId(*it) << " " << pos_it->second.size()
				<< endl;

		for (size_t i = 0; i < pos_it->second.size(); i++) {
			file << "    " << pos_it->second[i].contigId_ << ": "
					<< pos_it->second[i].start_ << " - "
					<< pos_it->second[i].end_ << endl;
		}
	}
}

template<class Graph>
class ConjugateDataPrinter: public DataPrinter<Graph> {
	typedef DataPrinter<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	ConjugateDataPrinter(Graph const& g, BaseIdTrackHandler<VertexId, EdgeId> const& int_ids) :
			base(g, int_ids) {
	}

	ConjugateDataPrinter(const GraphComponent<Graph>& graph_component,
			BaseIdTrackHandler<VertexId, EdgeId> const& int_ids) :
			base(GraphComponent<Graph>(graph_component, true), int_ids) {
	}

	template<class VertexIt>
	ConjugateDataPrinter(const Graph& g, VertexIt begin, VertexIt end,
			BaseIdTrackHandler<VertexId, EdgeId> const& int_ids) :
			base(GraphComponent<Graph>(g, begin, end, true), int_ids) {
	}

	std::string toPrint(VertexId v) const {
		stringstream ss;
		ss
				<< "Vertex "
				<< this->id_handler().ReturnIntId(v)
				<< " ~ "
				<< this->id_handler().ReturnIntId(
						this->component().g().conjugate(v)) << " .";
		return ss.str();
	}

	std::string toPrint(EdgeId e) const {
		stringstream ss;
		ss
				<< "Edge "
				<< this->id_handler().ReturnIntId(e)
				<< " : "
				<< this->id_handler().ReturnIntId(
						this->component().g().EdgeStart(e))
				<< " -> "
				<< this->id_handler().ReturnIntId(
						this->component().g().EdgeEnd(e))
				<< ", l = "
				<< this->component().g().length(e)
				<< " ~ "
				<< this->id_handler().ReturnIntId(
						this->component().g().conjugate(e)) << " .";
		return ss.str();
	}

};

template<class Graph>
class NonconjugateDataPrinter: public DataPrinter<Graph> {
	typedef DataPrinter<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	NonconjugateDataPrinter(Graph const& g,
			BaseIdTrackHandler<VertexId, EdgeId> const& int_ids) :
			base(g, int_ids) {
	}

	NonconjugateDataPrinter(const GraphComponent<Graph>& graph_component,
			BaseIdTrackHandler<VertexId, EdgeId> const& int_ids) :
			base(graph_component, int_ids) {
	}

	template<class VertexIt>
	NonconjugateDataPrinter(const Graph& g, VertexIt begin, VertexIt end,
			BaseIdTrackHandler<VertexId, EdgeId> const& int_ids) :
			base(GraphComponent<Graph>(g, begin, end), int_ids) {
	}

	std::string toPrint(VertexId v) const {
		stringstream ss;
		ss << "Vertex " << this->id_handler().ReturnIntId(v) << " .";
		return ss.str();
	}

	std::string toPrint(EdgeId e) const {
		stringstream ss;
		ss
				<< "Edge "
				<< this->id_handler().ReturnIntId(e)
				<< " : "
				<< this->id_handler().ReturnIntId(
						this->component().g().EdgeStart(e))
				<< " -> "
				<< this->id_handler().ReturnIntId(
						this->component().g().EdgeEnd(e)) << ", l = "
				<< this->component().g().length(e) << " .";
		return ss.str();
	}
};

template<class Graph>
struct PrinterTraits {
	typedef DataPrinter<Graph> Printer;
};

template<>
struct PrinterTraits<ConjugateDeBruijnGraph> {
	typedef ConjugateDataPrinter<ConjugateDeBruijnGraph> Printer;
};

template<>
struct PrinterTraits<NonconjugateDeBruijnGraph> {
	typedef NonconjugateDataPrinter<NonconjugateDeBruijnGraph> Printer;
};

template<class Graph>
class DataScanner {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	virtual void loadGraph(const string& file_name) = 0;
	void loadCoverage(const string& file_name);
	void loadPaired(const string& file_name,
			PairedInfoIndex<Graph>& paired_index);
	void loadPositions(const string& file_name,
			EdgesPositionHandler<Graph>& edge_pos);

private:
	Graph& g_;
//	int edge_count_;
	BaseIdTrackHandler<VertexId, EdgeId>& id_handler_;

protected:
	DataScanner(Graph &g, BaseIdTrackHandler<VertexId, EdgeId>& id_handler) :
			g_(g), id_handler_(id_handler) {
		INFO("Creating of scanner started");
//		edge_count_ = 0;
	}

	Graph& g() {
		return g_;
	}

	BaseIdTrackHandler<VertexId, EdgeId>& id_handler() {
		return id_handler_;
	}

public:
	virtual ~DataScanner() {

	}
};

template<class Graph>
class ConjugateDataScanner: public DataScanner<Graph> {
	typedef DataScanner<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	/*virtual*/
	void loadGraph(const string& file_name) {
		int flag;
		INFO(
				"Trying to read conjugate de bruijn  graph from " << file_name << ".grp");
		FILE* file = fopen((file_name + ".grp").c_str(), "r");
		VERIFY_MSG(file != NULL, "Couldn't find file " << (file_name + ".grp"));
		FILE* sequence_file = fopen((file_name + ".sqn").c_str(), "r");
		VERIFY(sequence_file != NULL);
		set<int> vertex_set;
		set<int> edge_set;
		INFO(
				"Reading conjugate de bruijn  graph from " << file_name << " started");
		size_t vertex_count;
		size_t edge_count;
		flag = fscanf(file, "%ld %ld \n", &vertex_count, &edge_count);
		VERIFY(flag == 2);
		for (size_t i = 0; i < vertex_count; i++) {
			size_t vertex_real_id, conjugate_id;
			flag = fscanf(file, "Vertex %ld ~ %ld .\n", &vertex_real_id,
					&conjugate_id);
			TRACE("Vertex "<<vertex_real_id<<" ~ "<<conjugate_id<<" .");
			VERIFY(flag == 2);

			if (vertex_set.find(vertex_real_id) == vertex_set.end()) {
				VertexId vid = this->g().AddVertex();
				VertexId conj_vid = this->g().conjugate(vid);

				this->id_handler().AddVertexIntId(vid, vertex_real_id);
				this->id_handler().AddVertexIntId(conj_vid, conjugate_id);
				vertex_set.insert(conjugate_id);
				TRACE(
						vid<<" ( "<< this->id_handler().ReturnVertexId(vertex_real_id) <<" )   "<< conj_vid << "( "<<this->id_handler().ReturnVertexId(conjugate_id)<<" )  added");
			}
		}

		char first_char = getc(sequence_file);
		VERIFY(!ferror(sequence_file));
		ungetc(first_char, sequence_file);
		bool fasta = (first_char == '>'); // if it's not fasta, then it's old .sqn


		if (!fasta) {
			size_t tmp_edge_count;
			flag = fscanf(sequence_file, "%ld", &tmp_edge_count);
			VERIFY(flag == 1);
			VERIFY(edge_count == tmp_edge_count);
		}

		const size_t longstring_size = 1000500; // TODO: O RLY magic constant? => Can't load edges >= 1Mbp
		char longstring[longstring_size];
		for (size_t i = 0; i < edge_count; i++) {
			size_t e_real_id, start_id, fin_id, length, conjugate_edge_id;
			flag = fscanf(file, "Edge %ld : %ld -> %ld, l = %ld ~ %ld .\n",
					&e_real_id, &start_id, &fin_id, &length,
					&conjugate_edge_id);
			VERIFY(flag == 5);
			VERIFY(length < longstring_size);
			if (fasta) {
				flag = fscanf(sequence_file, ">%ld\n%s\n", &e_real_id, longstring);
			}
			else {
				flag = fscanf(sequence_file, "%ld %s .", &e_real_id, longstring);
			}
			VERIFY(flag == 2);
			TRACE(
					"Edge "<<e_real_id<<" : "<<start_id<<" -> " << fin_id << " l = " << length << " ~ "<< conjugate_edge_id);
			if (edge_set.find(e_real_id) == edge_set.end()) {
				Sequence tmp(longstring);
				TRACE(
						start_id<<" "<< fin_id <<" "<< this->id_handler().ReturnVertexId(start_id)<<" "<< this->id_handler().ReturnVertexId(fin_id));
				EdgeId eid = this->g().AddEdge(
						this->id_handler().ReturnVertexId(start_id),
						this->id_handler().ReturnVertexId(fin_id), tmp);
				this->id_handler().AddEdgeIntId(eid, e_real_id);
				this->id_handler().AddEdgeIntId(this->g().conjugate(eid),
						conjugate_edge_id);
				edge_set.insert(conjugate_edge_id);
			}
		}
		fclose(file);
		fclose(sequence_file);
	}
public:
	ConjugateDataScanner(Graph& g, BaseIdTrackHandler<VertexId, EdgeId>& id_handler) :
			base(g, id_handler) {
	}
};

template<class Graph>
class NonconjugateDataScanner: public DataScanner<Graph> {
	typedef DataScanner<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	/*virtual*/
	void loadGraph(const string& file_name) {
		int flag;
		FILE* file = fopen((file_name + ".grp").c_str(), "r");
		if (file == NULL
)
						WARN("File "<<(file_name + ".grp")<<" not found");
		VERIFY(file != NULL);
		FILE* sequence_file = fopen((file_name + ".sqn").c_str(), "r");
		VERIFY(sequence_file != NULL);

		INFO(
				"Reading NON conjugate de bruujn graph from " << file_name << " started");
		size_t vertex_count;
		size_t edge_count;
		flag = fscanf(file, "%ld %ld \n", &vertex_count, &edge_count);
		VERIFY(flag == 2);
		for (size_t i = 0; i < vertex_count; i++) {
			size_t vertex_real_id;
			flag = fscanf(file, "Vertex %ld", &vertex_real_id);
			VERIFY(flag == 1);
			char c = 'a';
			while (c != '.') {
				flag = fscanf(file, "%c", &c);
				VERIFY(flag == 1);
			}
			flag = fscanf(file, "\n");
			VERIFY(flag == 0);
			VertexId vid = this->g().AddVertex();
			this->id_handler().AddVertexIntId(vid, vertex_real_id);
			TRACE(vid);
		}
		size_t tmp_edge_count;
		flag = fscanf(sequence_file, "%ld", &tmp_edge_count);
		VERIFY(flag == 1);
		VERIFY(edge_count == tmp_edge_count);
		char longstring[1000500];
		for (size_t i = 0; i < edge_count; i++) {
			int e_real_id, start_id, fin_id, length;
			flag = fscanf(file, "Edge %d : %d -> %d, l = %d", &e_real_id,
					&start_id, &fin_id, &length);
			VERIFY(flag == 4);
			flag = fscanf(sequence_file, "%d %s .", &e_real_id, longstring);
			VERIFY(flag == 2);
			//does'nt matter, whether it was conjugate or not.
			char c = 'a';
			while (c != '.') {
				flag = fscanf(file, "%c", &c);
				VERIFY(flag == 1);
			}
			flag = fscanf(file, "\n");
			VERIFY(flag == 0);
			Sequence tmp(longstring);
			TRACE(
					start_id<<" "<< fin_id <<" "<< this->id_handler().ReturnVertexId(start_id)<<" "<< this->id_handler().ReturnVertexId(fin_id));
			EdgeId eid = this->g().AddEdge(
					this->id_handler().ReturnVertexId(start_id),
					this->id_handler().ReturnVertexId(fin_id), tmp);
			this->id_handler().AddEdgeIntId(eid, e_real_id);

		}
		fclose(file);
		fclose(sequence_file);
	}

	NonconjugateDataScanner(Graph &g, BaseIdTrackHandler<VertexId, EdgeId>& id_handler) :
			base(g, id_handler) {
	}
};

template<class Graph>
void DataScanner<Graph>::loadCoverage(const string& file_name) {
	int read_count;
	FILE* file = fopen((file_name + ".cvr").c_str(), "r");
	VERIFY(file != NULL);
	INFO("Reading coverage from " << file_name << " started");
	int edge_count;
	read_count = fscanf(file, "%d \n", &edge_count);
	VERIFY(read_count == 1);
//	VERIFY(edge_count == edge_count_);
	for (int i = 0; i < edge_count; i++) {
		int edge_real_id;
		double edge_coverage;
		read_count = fscanf(file, "%d %lf .\n", &edge_real_id, &edge_coverage);
		VERIFY(read_count == 2);
		TRACE(edge_real_id<< " "<<edge_coverage <<" . ");
		EdgeId eid = id_handler_.ReturnEdgeId(edge_real_id);
		TRACE("EdgeId "<<eid);
		g_.coverage_index().SetCoverage(eid, edge_coverage * g_.length(eid));
	}
	fclose(file);
}

template<class Graph>
void DataScanner<Graph>::loadPaired(const string& file_name,
		PairedInfoIndex<Graph>& paired_index) {
	typedef typename Graph::EdgeId EdgeId;
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
		TRACE(
				first_real_id<< " " << second_real_id << " " << d << " " << w << " " << v);
		if (id_handler_.ReturnEdgeId(first_real_id) == EdgeId(NULL) || id_handler_.ReturnEdgeId(second_real_id) == EdgeId(NULL))
			continue;
		TRACE(
				id_handler_.ReturnEdgeId(first_real_id)<<" "<< id_handler_.ReturnEdgeId(second_real_id)<<" "<< d<<" "<< w);
		PairInfo<typename Graph::EdgeId> p_info(
				id_handler_.ReturnEdgeId(first_real_id),
				id_handler_.ReturnEdgeId(second_real_id), d, w, v);
		paired_index.AddPairInfo(p_info, false);
		PairInfo<typename Graph::EdgeId> p_info(
				id_handler_.ReturnEdgeId(second_real_id),
				id_handler_.ReturnEdgeId(first_real_id), -d, w, v);
		paired_index.AddPairInfo(p_info, false);
	}DEBUG("PII SIZE " << paired_index.size());
	fclose(file);
}

template<class Graph>
void DataScanner<Graph>::loadPositions(const string& file_name,
		EdgesPositionHandler<Graph>& edge_pos) {
	int read_count;
	FILE* file = fopen((file_name + ".pos").c_str(), "r");
	VERIFY(file != NULL);
	DEBUG("Reading edges positions, " << file_name <<" started");
	VERIFY(file != NULL);
	int pos_count;
	read_count = fscanf(file, "%d\n", &pos_count);
	VERIFY(read_count == 1);
	for (int i = 0; i < pos_count; i++) {
		int edge_real_id, pos_info_count;
		char contigId[500];
		read_count = fscanf(file, "%d %d\n", &edge_real_id, &pos_info_count);
		VERIFY(read_count == 2);
//		INFO(  edge_real_id);
		for (int j = 0; j < pos_info_count; j++) {
			int start_pos, end_pos;
			read_count = fscanf(file, "%s %d - %d \n", contigId, &start_pos,
					&end_pos);
//			INFO (contigId<<" "<< start_pos<<" "<<end_pos);
			VERIFY(read_count == 3);
			EdgeId eid = id_handler_.ReturnEdgeId(edge_real_id);
			edge_pos.AddEdgePosition(eid, start_pos, end_pos, string(contigId));
		}
	}
	fclose(file);
}

template<class Graph>
struct ScannerTraits {
	typedef DataScanner<Graph> Scanner;
};

template<>
struct ScannerTraits<ConjugateDeBruijnGraph> {
	typedef ConjugateDataScanner<ConjugateDeBruijnGraph> Scanner;
};

template<>
struct ScannerTraits<NonconjugateDeBruijnGraph> {
	typedef NonconjugateDataScanner<NonconjugateDeBruijnGraph> Scanner;
};

//helper methods
// todo think how to organize them in the most natural way

template<class Graph>
void PrintBasicGraph(const string& file_name, DataPrinter<Graph>& printer) {
	printer.saveGraph(file_name);
	printer.saveEdgeSequences(file_name);
	printer.saveCoverage(file_name);
}

template<class graph_pack>
void PrintGraphPack(const string& file_name,
		DataPrinter<typename graph_pack::graph_t>& printer,
		const graph_pack& gp) {
	PrintBasicGraph(file_name, printer);
	printer.savePaired(file_name + "_et", gp.etalon_paired_index);
	printer.savePositions(file_name, gp.edge_pos);
	SaveKmerMapper(file_name, gp.kmer_mapper);
}

template<class graph_pack>
void PrintGraphPack(const string& file_name, const graph_pack& gp) {
	typename PrinterTraits<typename graph_pack::graph_t>::Printer printer(gp.g,
			gp.int_ids);
	PrintGraphPack(file_name, printer, gp);
}

template<class Graph>
void PrintPairedIndex(const string& file_name, DataPrinter<Graph>& printer,
		const PairedInfoIndex<Graph>& paired_index) {
	printer.savePaired(file_name, paired_index);
}

template<class Graph>
void PrintClusteredIndex(const string& file_name, DataPrinter<Graph>& printer,
		const PairedInfoIndex<Graph>& clustered_index) {
	PrintPairedIndex(file_name + "_cl", printer, clustered_index);
}

template<class graph_pack>
void PrintWithPairedIndex(const string& file_name,
		DataPrinter<typename graph_pack::graph_t>& printer,
		const graph_pack& gp,
		const PairedInfoIndex<typename graph_pack::graph_t>& paired_index,
		bool clustered_index = false) {
	PrintGraphPack(file_name, printer, gp);
	if (!clustered_index) {
		PrintPairedIndex(file_name, printer, paired_index);
	} else {
		PrintClusteredIndex(file_name, printer, paired_index);
	}
}

template<class graph_pack>
void PrintWithClusteredIndex(const string& file_name,
		DataPrinter<typename graph_pack::graph_t>& printer,
		const graph_pack& gp,
		const PairedInfoIndex<typename graph_pack::graph_t>& paired_index) {
	PrintWithPairedIndex(file_name, printer, gp, paired_index, true);
}

template<class graph_pack, class VertexIt>
void PrintAll(const string& file_name, const graph_pack& gp, VertexIt begin,
		VertexIt end,
		const PairedInfoIndex<typename graph_pack::graph_t>& paired_index,
		const PairedInfoIndex<typename graph_pack::graph_t>& clustered_index) {
	typename PrinterTraits<typename graph_pack::graph_t>::Printer printer(gp.g,
			begin, end, gp.int_ids);
	PrintGraphPack(file_name, printer, gp);
	PrintPairedIndex(file_name, printer, paired_index);
	PrintClusteredIndex(file_name, printer, clustered_index);
}

template<class graph_pack>
void PrintAll(const string& file_name, const graph_pack& gp,
		const PairedInfoIndex<typename graph_pack::graph_t>& paired_index,
		const PairedInfoIndex<typename graph_pack::graph_t>& clustered_index) {
	PrintAll(file_name, gp, gp.g.begin(), gp.g.end(), paired_index,
			clustered_index);
}

template<class graph_pack, class VertexIt>
void PrintWithPairedIndex(const string& file_name, const graph_pack& gp,
		VertexIt begin, VertexIt end,
		const PairedInfoIndex<typename graph_pack::graph_t>& paired_index,
		bool clustered_index = false) {
	typename PrinterTraits<typename graph_pack::graph_t>::Printer printer(gp.g,
			begin, end, gp.int_ids);
	PrintWithPairedIndex(file_name, printer, gp, paired_index, clustered_index);
}

template<class graph_pack, class VertexIt>
void PrintWithClusteredIndex(const string& file_name, const graph_pack& gp,
		VertexIt begin, VertexIt end,
		const PairedInfoIndex<typename graph_pack::graph_t>& clustered_index) {
	typename PrinterTraits<typename graph_pack::graph_t>::Printer printer(gp.g,
			begin, end, gp.int_ids);
	PrintWithPairedIndex(file_name, printer, gp, clustered_index, true);
}

template<class graph_pack>
void PrintWithPairedIndex(const string& file_name, const graph_pack& gp,
		const PairedInfoIndex<typename graph_pack::graph_t>& paired_index,
		bool clustered_index = false) {
	PrintWithPairedIndex(file_name, gp, gp.g.begin(), gp.g.end(), paired_index,
			clustered_index);
}

template<class graph_pack, class VertexIt>
void PrinGraphPack(const string& file_name, const graph_pack& gp,
		VertexIt begin, VertexIt end) {
	typename PrinterTraits<typename graph_pack::graph_t>::Printer printer(gp.g,
			begin, end, gp.int_ids);
	PrintGraphPack(file_name, printer, gp);
}

template<class graph_pack>
void PrintWithClusteredIndex(const string& file_name, const graph_pack& gp,
		const PairedInfoIndex<typename graph_pack::graph_t>& clustered_index) {
	PrintWithPairedIndex(file_name, gp, clustered_index, true);
}

template<class Graph>
void ScanBasicGraph(const string& file_name, DataScanner<Graph>& scanner) {
	scanner.loadGraph(file_name);
	scanner.loadCoverage(file_name);
}

template<class graph_pack>
void ScanGraphPack(const string& file_name,
		DataScanner<typename graph_pack::graph_t>& scanner, graph_pack& gp) {
	ScanBasicGraph(file_name, scanner);
	scanner.loadPaired(file_name + "_et", gp.etalon_paired_index);
	scanner.loadPositions(file_name, gp.edge_pos);
	LoadKmerMapper(file_name, gp.kmer_mapper);
}

template<class Graph>
void ScanPairedIndex(const string& file_name, DataScanner<Graph>& scanner,
		PairedInfoIndex<Graph>& paired_index) {
	scanner.loadPaired(file_name, paired_index);
}

template<class Graph>
void ScanClusteredIndex(const string& file_name, DataScanner<Graph>& scanner,
		PairedInfoIndex<Graph>& clustered_index) {
	scanner.loadPaired(file_name + "_cl", clustered_index);
}

template<class graph_pack>
void ScanWithPairedIndex(const string& file_name,
		DataScanner<typename graph_pack::graph_t>& scanner, graph_pack& gp,
		PairedInfoIndex<typename graph_pack::graph_t>& paired_index,
		bool clustered_index = false) {
	ScanGraphPack(file_name, scanner, gp);
	if (!clustered_index) {
		ScanPairedIndex(file_name, scanner, paired_index);
	} else {
		ScanClusteredIndex(file_name, scanner, paired_index);
	}
}

template<class graph_pack>
void ScanWithClusteredIndex(const string& file_name,
		DataScanner<typename graph_pack::graph_t>& scanner, graph_pack& gp,
		PairedInfoIndex<typename graph_pack::graph_t>& paired_index) {
	ScanWithPairedIndex(file_name, scanner, gp, paired_index, true);
}

template<class graph_pack>
void ScanWithPairedIndex(const string& file_name, graph_pack& gp,
		PairedInfoIndex<typename graph_pack::graph_t>& paired_index,
		bool clustered_index = false) {
	typename ScannerTraits<typename graph_pack::graph_t>::Scanner scanner(gp.g,
			gp.int_ids);
	ScanWithPairedIndex(file_name, scanner, gp, paired_index, clustered_index);
}

template<class Graph>
void ScanBasicGraph(const string& file_name, Graph& g,
		BaseIdTrackHandler<VertexId, EdgeId>& int_ids) {
	typename ScannerTraits<Graph>::Scanner scanner(g, int_ids);
	ScanBasicGraph(file_name, scanner);
}

template<class graph_pack>
void ScanGraphPack(const string& file_name, graph_pack& gp) {
	typename ScannerTraits<typename graph_pack::graph_t>::Scanner scanner(gp.g,
			gp.int_ids);
	ScanGraphPack(file_name, scanner, gp);
}

template<class graph_pack>
void ScanAll(const string& file_name, graph_pack& gp,
		PairedInfoIndex<typename graph_pack::graph_t>& paired_index,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index) {
	typename ScannerTraits<typename graph_pack::graph_t>::Scanner scanner(gp.g,
			gp.int_ids);
	ScanGraphPack(file_name, scanner, gp);
	ScanPairedIndex(file_name, scanner, paired_index);
	ScanClusteredIndex(file_name, scanner, clustered_index);
}

template<class graph_pack>
void ScanWithClusteredIndex(const string& file_name, graph_pack& gp,
		PairedInfoIndex<typename graph_pack::graph_t>& clustered_index) {
	ScanWithPairedIndex(file_name, gp, clustered_index, true);
}

}
