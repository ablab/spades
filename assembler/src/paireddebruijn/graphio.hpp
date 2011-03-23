#ifndef IOPROCEDURES_HPP_
#define IOPROCEDURES_HPP_
#include "common.hpp"
#include "pairedGraph.hpp"
using namespace paired_assembler;

void outputLongEdges(longEdgesMap &longEdges);
void outputLongEdges(longEdgesMap &longEdges, PairedGraph &graph);
void outputLongEdgesThroughGenome(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount);

void codeRead(char *read, char *code);

inline bool nextReadPair(char * &read1, char * &read2) {
	return (scanf("%s %s", read1, read2) == 2);
}

ll extractMer(char *read, int shift, int length);

string decompress(ll a, int l);

class DataPrinter {
	FILE *f_;
public:
	DataPrinter(char *fileName);
	void output(int a);
	void output(long long a);
	void output(Edge *edge);
	void output(Sequence *sequence);
	void output(VertexPrototype *v);
	void outputLongEdgesMap(longEdgesMap &map);
	void outputIntArray(int *array, int length);
	void outputIntArray(int *array, int length, int width);
	template<typename keyType, typename valueType>
	void output(map<keyType, valueType> m);
	template<typename valueType>
	void output(vector<valueType> v);
	void close();
};

class DataReader {
	FILE *f_;
public:
	DataReader(char *fileName);
	void read(int &a);
	void read(long long &a);
	void read(Edge * &edge);
	void read(Sequence * &sequence);
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

void save(char *fileName, PairedGraph &g, longEdgesMap &longEdges,
		int &VertexCount, int EdgeId);

void load(char *fileName, PairedGraph &g, longEdgesMap &longEdges,
		int &VertexCount, int EdgeId);

#endif /* IOPROCEDURES_HPP_ */
