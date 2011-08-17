

//#include "graphio.hpp"
////#include "graphSimplification.hpp"
////#include "graphVisualizer.hpp"
//#include <stdio.h>
//#include "common.hpp"
//#include "iostream"
//#include "fstream"
//
//DECL_MODULE_LOGGER("graphio");
//using namespace omni;
//
//DataReader::DataReader(char *fileName) {
//	f_ = fopen(fileName, "r");
//	DEBUG("DataReader " << fileName <<" created");
//	assert(f_ != NULL);
//}
//
//DataReader::DataReader(const char *fileName) {
//	f_ = fopen(fileName, "r");
//	DEBUG("DataReader " << fileName <<" created");
//	assert(f_ != NULL);
//}
//
//DataPrinter::DataPrinter(char *fileName) {
//	f_ = fopen(fileName, "w");
//	assert(f_ != NULL);
//}
//
//DataPrinter::DataPrinter(const char *fileName) {
//	f_ = fopen(fileName, "w");
//	assert(f_ != NULL);
//}
//
//void DataReader::close() {
//	fclose(f_);
//}
//
//void DataPrinter::close() {
//	fclose(f_);
//}
//
//int DataReader::read(int &a) {
//	return fscanf(f_, "%d\n", &a);
//}
//
//void DataPrinter::output(int a) {
//	fprintf(f_, "%d\n", a);
//}
//
//int DataReader::read(long long &a) {
//	return fscanf(f_, "%lld\n", &a);
//}
//
//void DataPrinter::output(long long a) {
//	fprintf(f_, "%lld\n", a);
//}
//
//int DataReader::read(Sequence * &sequence) {
//	int length;
//	if (!read(length)) {
//		return 0;
//	} else {
//		if (length == 0) {
//			fscanf(f_, "\n");
//			sequence = new Sequence("");
//		} else {
//			char *s = new char[length + 1];
//			fscanf(f_, "%s\n", s);
//			sequence = new Sequence(s);
//		}
//		return 1;
//	}
//}
//
//void DataPrinter::output(Sequence *sequence) {
//	output((int) sequence->size());
//	fprintf(f_, "%s\n", sequence->str().c_str());
//}
//
//void DataReader::read(VertexPrototype * &v) {
//	int id;
//	Sequence *lower;
//	bool b;
//	read(id);
//	read(lower);
//	int tmpInt;
//	read(tmpInt);
//	b = tmpInt;
//	v = new VertexPrototype(lower, id);
//	v->used = b;
//}
//
//void DataPrinter::output(VertexPrototype *v) {
//	output(v->VertexId);
//	output(v->lower);
//	output(v->used);
//}
//
//int DataReader::read(Edge * &edge) {
//	int from, to, len, id, cov;
//	Sequence *up, *low;
//	int read_res = 0;
//	read_res += read(id);
//	if (id == -1) return 0;
//	DEBUG(id);
//	read_res += read(from);
//	read_res += read(to);
//	read_res += read(len);
//	read_res += read(cov);
//	read_res += read(up);
//	read_res += read(low);
//	edge = new Edge(up, low, from, to, len, id, cov);
//	DEBUG("edge"<<id << "loaded");
//	return (read_res/7);
//}
//
//void DataPrinter::output(Edge *edge) {
//	output(edge->EdgeId);
//	output(edge->FromVertex);
//	output(edge->ToVertex);
//	output(edge->length);
//	output(edge->coverage);
//	output(edge->upper);
//	output(edge->lower);
//}
//
//void DataPrinter::outputLongEdgesMap(longEdgesMap &edges) {
//	INFO("Saving long edges");
//	DEBUG(edges.size());
//	int size = edges.size();
//	for (longEdgesMap::iterator it = edges.begin(); it != edges.end(); ++it) {
//		if (it->first == it->second->EdgeId) {
//			output(it->first);
//			output(it->second);
//			DEBUG("Edge outputed" );
//		}
//	}
//	DEBUG("Normal edges outputed");
//	Sequence *emptySequence = new Sequence("");
//	Edge *emptyEdge = new Edge(emptySequence, emptySequence, 0, 0, 0, 0);
//	for (longEdgesMap::iterator it = edges.begin(); it != edges.end(); ++it) {
//		if (it->first != it->second->EdgeId) {
//			output(it->first);
//			emptyEdge->EdgeId = it->second->EdgeId;
//			output(emptyEdge);
//		}
//	}
//	DEBUG("fakeEdges outputed");
//	delete emptyEdge;
//	output(-1);
//	output(size);
//}
//
//void DataReader::readLongEdgesMap(longEdgesMap &edges) {
//	int size = 0;
//	int id;
//	while (1) {
//		assert(read(id));
//		Edge *edge;
//		DEBUG(id);
//		if( id == -1 || !(read(edge)))
//			break;
//		else {
//			size++;
//			if (id == edge->EdgeId) {
//				edges.insert(make_pair(id, edge));
//			} else {
//				edges.insert(make_pair(id, edges[edge->EdgeId]));
//				delete edge;
//			}
//		}
//	}
//	read(id);
//	DEBUG(id);
//	assert(size == id);
//}
//
//void DataReader::readIntArray(int *array, int length) {
//	for (int i = 0; i < length; i++) {
//		fscanf(f_, "%d ", array + i);
//	}
//	fscanf(f_, "\n");
//}
//
//void DataPrinter::outputIntArray(int *array, int length) {
//	for (int i = 0; i < length; i++) {
//		fprintf(f_, "%d ", array[i]);
//	}
//	fprintf(f_, "\n");
//}
//
//void DataReader::readIntArray(int *array, int length, int width) {
//	int cur = 0;
//	for (int i = 0; i < length; i++) {
//		for (int j = 0; j < width; j++) {
//			fscanf(f_, "%d ", array + cur);
//			cur++;
//		}
//		fscanf(f_, "\n");
//	}
//	fscanf(f_, "\n");
//}
//
//void DataPrinter::outputIntArray(int *array, int length, int width) {
//	int cur = 0;
//	for (int i = 0; i < length; i++) {
//		for (int j = 0; j < width; j++) {
//			fprintf(f_, "%d ", array[cur]);
//			cur++;
//		}
//		fprintf(f_, "\n");
//	}
//	fprintf(f_, "\n");
//}
//
//void save(DataPrinter dp, Edge *e){
//	dp.output(e->EdgeId);
//	dp.output(e);
////	dp.close();
//}
//void save(char *fileName, PairedGraph &g, longEdgesMap &longEdges,
//		int &VertexCount, int EdgeId) {
//	DataPrinter dp(fileName);
//	dp.outputLongEdgesMap(longEdges);
//	dp.output(VertexCount);
//	dp.output(EdgeId);
//	//TODO: FIX!!!
//	//	dp.outputIntArray(g.inD, MAX_VERT_NUMBER);
//	//	dp.outputIntArray(g.outD, MAX_VERT_NUMBER);
//	//	dp.outputIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	//	dp.outputIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dp.close();
//}
//void load(char *fileName, PairedGraph &g, longEdgesMap &longEdges,
//		int &VertexCount, int &EdgeId) {
//	DataReader dr(fileName);
//	dr.readLongEdgesMap(longEdges);
//	dr.read(VertexCount);
//	dr.read(EdgeId);
//	//TODO: fix;
//	//	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
//	//	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
//	//	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	//	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dr.close();
//}
//
//void save(char *fileName, PairedGraph &g) {
//	INFO("Saving graph");
//	DataPrinter dp(fileName);
//	dp.output(g.VertexCount);
//	dp.output(g.EdgeId);
//	dp.outputLongEdgesMap(g.longEdges);
//	//TODO: FIX!!!
//	//	dp.outputIntArray(g.inD, MAX_VERT_NUMBER);
//	//	dp.outputIntArray(g.outD, MAX_VERT_NUMBER);
//	//	dp.outputIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	//	dp.outputIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dp.close();
//}
//void save(string fileName, PairedGraph &g) {
//	INFO("Saving graph");
//	DataPrinter dp(fileName.c_str());
//	dp.outputLongEdgesMap(g.longEdges);
//	dp.output(g.VertexCount);
//	dp.output(g.EdgeId);
//	//TODO: FIX!!!
//	//	dp.outputIntArray(g.inD, MAX_VERT_NUMBER);
//	//	dp.outputIntArray(g.outD, MAX_VERT_NUMBER);
//	//	dp.outputIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	//	dp.outputIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dp.close();
//}
//void load(char *fileName, PairedGraph &g) {
//	INFO("Loading graph");
//	DataReader dr(fileName);
//	dr.readLongEdgesMap(g.longEdges);
//	dr.read(g.VertexCount);
//	dr.read(g.EdgeId);
//	//TODO: fix;
//	//	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
//	//	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
//	//	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	//	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dr.close();
//}
//void load(string fileName, PairedGraph &g) {
//	INFO("Loading graph");
//	DataReader dr(fileName.c_str());
//	dr.readLongEdgesMap(g.longEdges);
//	dr.read(g.VertexCount);
//	dr.read(g.EdgeId);
//	//TODO: fix;
//	//	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
//	//	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
//	//	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	//	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dr.close();
//}
//
//void outputVertexKmers(edgesMap &edges){
//
//	FILE *fkmers = fopen((folder+"kmers.txt").c_str(), "w");
//	for (edgesMap::iterator iter = edges.begin(); iter != edges.end();++iter) {
//		ll kmer = iter->fi;
//		fprintf(fkmers,"%lld %s\n", kmer, decompress(kmer, k));
//	}
//	fclose(fkmers);
//
//}
