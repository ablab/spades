#ifndef IOPROCEDURES_HPP_
#define IOPROCEDURES_HPP_
#include "common.hpp"
#include "pairedGraph.hpp"
using namespace paired_assembler;

void outputLongEdges(longEdgesMap &longEdges);
void outputLongEdges(longEdgesMap &longEdges, string fileName);
void outputLongEdges(longEdgesMap &longEdges, PairedGraph &graph);
void outputLongEdges(longEdgesMap &longEdges, PairedGraph &graph, string fileName);
void outputLongEdgesThroughGenome(PairedGraph &graph);
void outputLongEdgesThroughGenome(PairedGraph &graph, string fileName);

void codeRead(char *read, char *code);

void outputVertexKmers(PairedGraph &graph);

inline bool nextReadPair(FILE* f, char * &read1, char * &read2) {
	//if (!fictiveSecondReads
   return (fscanf(f, "%s %s", read1, read2) == 2);
/*	else {
		if (fscanf(f, "%s %s", read1, read2) == 2){
			forn(i,strlen(read2)) read2[i]='A';
			return true;
		}
		return false;
	}
	*/
}

Sequence ReadGenome(istream &is);

Sequence ReadGenomeFromFile(const string &fileName);


ll extractMer(char *read, int shift, int length);

string decompress(ll a, int l);

class DataPrinter {
	FILE *f_;
public:
	DataPrinter(char *fileName);
	DataPrinter(const char *fileName);
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

private:
	DECL_LOGGER("DataPrinter")
};

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

void save(char *fileName, PairedGraph &g, longEdgesMap &longEdges,
		int &VertexCount, int EdgeId);
void save(DataPrinter dp, Edge *e);
void load(char *fileName, PairedGraph &g, longEdgesMap &longEdges,
		int &VertexCount, int EdgeId);

void save(char *fileName, PairedGraph &g);
void load(char *fileName, PairedGraph &g);

void save(string fileName, PairedGraph &g);
void load(string fileName, PairedGraph &g);

#endif /* IOPROCEDURES_HPP_ */
