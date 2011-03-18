#include "graphio.hpp"
#include "graphVisualizer.hpp"
#include "common.hpp"
#include "pairedGraph.hpp"

using namespace paired_assembler;

inline int codeNucleotide(char a) {
	if (a == 'A')
		return 0;
	else if (a == 'C')
		return 1;
	else if (a == 'G')
		return 2;
	else if (a == 'T')
		return 3;
	else {
		std::cerr << "oops!";
		return -1;
	}
}

void codeRead(char *read, char *code) {
	for (int i = 0; i < readLength; i++) {
		code[i] = codeNucleotide(read[i]);
	}
}

ll extractMer(char *read, int shift, int length) {
	ll res = 0;
	for (int i = 0; i < length; i++) {
		res = res << 2;
		res += read[shift + i];
	}
	return res;
}

string decompress(ll a, int l) {

	string res = "";
	res.reserve(l);
	forn(i,l)
		res += " ";
	forn(i, l) {
		res[l - i - 1] = nucl((a & 3));
		a >>= 2;
	}
	return res;
}

void outputLongEdges(longEdgesMap &longEdges) {
	char Buffer[100];
	gvis::GraphPrinter<int> g("Paired_ext");
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			sprintf(Buffer, "%i (%i)", it->first, it->second->length);
			//		else sprintf(Buffer,"%i (%i) FAKE now it is %d",it->first, it->second->length,it->second->EdgeId);

			g.addEdge(it->second->FromVertex, it->second->ToVertex, Buffer);
			cerr << it->first << " (" << it->second->length << "):" << endl;
			cerr << it->second->upper->str() << endl;
			cerr << it->second->lower->str() << endl;
		}
	}
	g.output();
}

DataPrinter::DataPrinter(char *fileName) {
	f_ = fopen(fileName, "w");
}

void DataPrinter::close() {
	fclose(f_);
}

void DataPrinter::outputInt(int a) {
	fprintf(f_, "%d\n", a);
}

void DataPrinter::outputSequence(Sequence *sequence) {
	outputInt(sequence->size());
	fprintf(f_, "%s\n", sequence->str().c_str());
}

void DataPrinter::outputEdge(Edge *edge) {
	outputInt(edge->EdgeId);
	outputInt(edge->FromVertex);
	outputInt(edge->ToVertex);
	outputInt(edge->length);
	outputSequence(edge->upper);
	outputSequence(edge->lower);
}

void DataPrinter::outputLongEdgesMap(longEdgesMap &edges) {
	outputInt(edges.size());
	for(longEdgesMap::iterator it = edges.begin(); it != edges.end(); ++it) {
		outputInt(it->first);
		outputEdge(it->second);
	}
}

void DataPrinter::outputIntArray(int *array, int length) {
	for (int i = 0; i < length; i++) {
		fprintf(f_, "%d ", array[i]);
	}
	fprintf(f_, "\n");
}

void DataPrinter::outputIntArray(int *array, int length, int width) {
	int cur = 0;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			fprintf(f_, "%d ", array[cur]);
			cur++;
		}
		fprintf(f_, "\n");
	}
	fprintf(f_, "\n");
}

DataReader::DataReader(char *fileName) {
	f_ = fopen(fileName, "r");
}

void DataReader::close() {
	fclose(f_);
}

void DataReader::readInt(int &a) {
	fscanf(f_, "%d\n", &a);
}

void DataReader::readSequence(Sequence * &sequence) {
	int length;
	readInt(length);
	char *s = new char[length + 1];
	fscanf(f_, "%s\n", s);
	sequence = new Sequence(s);
}

void DataReader::readEdge(Edge * &edge){
	int from, to, len, id;
	Sequence *up, *low;
	readInt(id);
	readInt(from);
	readInt(to);
	readInt(len);
	readSequence(up);
	readSequence(low);
	edge = new Edge(up, low, from, to, len, id);
}

void DataReader::readLongEdgesMap(longEdgesMap &edges) {
	int size;
	readInt(size);
	for(int i = 0; i < size; i++) {
		int id;
		readInt(id);
		Edge *edge;
		readEdge(edge);
		edges.insert(make_pair(id, edge));
	}
}

void DataReader::readIntArray(int *array, int length) {
	for (int i = 0; i < length; i++) {
		fscanf(f_, "%d ", array + i);
	}
	fscanf(f_, "\n");
}

void DataReader::readIntArray(int *array, int length, int width) {
	int cur = 0;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			fscanf(f_, "%d ", array + cur);
			cur++;
		}
		fscanf(f_, "\n");
	}
	fscanf(f_, "\n");
}

void save(char *fileName, PairedGraph &g, longEdgesMap &longEdges, int &VertexCount,
		int EdgeId) {
	DataPrinter dp(fileName);
	dp.outputInt(VertexCount);
	dp.outputInt(EdgeId);
	dp.outputLongEdgesMap(longEdges);
	dp.outputIntArray(g.inD, MAX_VERT_NUMBER);
	dp.outputIntArray(g.outD, MAX_VERT_NUMBER);
	dp.outputIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dp.outputIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dp.close();
}
void load(char *fileName, PairedGraph &g, longEdgesMap &longEdges, int &VertexCount,
		int EdgeId) {
	DataReader dr(fileName);
	dr.readInt(VertexCount);
	dr.readInt(EdgeId);
	dr.readLongEdgesMap(longEdges);
	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dr.close();
}
