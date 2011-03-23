#include "graphio.hpp"
#include "graphVisualizer.hpp"
#include <stdio.h>
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
			if (it->second->length < 500) {
				cerr << it->second->upper->str() << endl;
				cerr << it->second->lower->str() << endl;
			}
		}
	}
	g.output();
}

void outputLongEdgesThroughGenome(longEdgesMap &longEdges, PairedGraph &graph,
		int &VertexCount) {
	char Buffer[100];
	char* Genome;
	int GenLength;
	int GenPos;
	int i, t;

	int EdgeNum = 0;

	int bigShift = insertLength + readLength;
	assert(k==l);

	cerr << "Graph output through genome" << endl;
	gvis::GraphPrinter<int> g("Paired_ext");

	FILE *infile;
	Genome = (char *) malloc(6000000);
	if ((infile = fopen("data/MG1655-K12_cut.fasta", "r")) == NULL) {
		cerr << "No such file" << endl;
		return;
	}
	i = 0;
	Genome[i] = fgetc(infile);
	t = 0;
	while (t < 5) {
		if (Genome[i] > 20) {
			t = 0;
			i++;
		} else
			t++;
		Genome[i] = fgetc(infile);
	}
	GenLength = i;
	GenPos = 0;
	fclose(infile);
	cerr << "Try to process" << endl;

	int CurVert = 0;
	while ((graph.degrees[CurVert][0]!=0)||(graph.degrees[CurVert][1]!=1)) CurVert++;
	cerr<<"Start vertex "<<CurVert<<endl;
	while (graph.degrees[CurVert][1]!=0){
		bool NoEdge = true;
		cerr<<"Try to found next edge"<<endl;
		forn(v,graph.degrees[CurVert][1])
		{
			int edgeId = edgeRealId(graph.edgeIds[CurVert][v][OUT_EDGE], longEdges);
			cerr << "possible edge" << edgeId << endl;
			bool goodEdge = true;
			int h = 0;
			while (goodEdge && (h < longEdges[edgeId]->upper->size())) {
				//cerr<<" "<<(*(longEdges[edgeId]->upper))[h]<<" =?= "<<Genome[GenPos+h]<<endl;
				if (nucl((*(longEdges[edgeId]->upper))[h])
						!= Genome[GenPos + h])
					goodEdge = false;
				h++;
			}
			h = 0;
			while (goodEdge && (h < longEdges[edgeId]->lower->size())) {
				if (nucl((*(longEdges[edgeId]->lower))[h]) != Genome[GenPos + h
						+ bigShift])
					goodEdge = false;
				h++;
			}

			if (goodEdge) {
				cerr << "Edge found" << endl;
				EdgeNum++;
				sprintf(Buffer, "%i: %i (%i)", EdgeNum, edgeId,
						longEdges[edgeId]->length);
				g.addEdge(longEdges[edgeId]->FromVertex,
						longEdges[edgeId]->ToVertex, Buffer);
				cerr << edgeId << " (" << longEdges[edgeId]->length << "):"
						<< endl;
				if (longEdges[edgeId]->length < 500) {
					cerr << longEdges[edgeId]->upper->str() << endl;
					cerr << longEdges[edgeId]->lower->str() << endl;
				}
				CurVert = longEdges[edgeId]->ToVertex;
				GenPos += longEdges[edgeId]->length;
				NoEdge = false;
				break;
			}
		}
		if (NoEdge)
			break;
	}
	g.output();
}

DataReader::DataReader(char *fileName) {
	f_ = fopen(fileName, "r");
}

DataPrinter::DataPrinter(char *fileName) {
	f_ = fopen(fileName, "w");
}

void DataReader::close() {
	fclose(f_);
}

void DataPrinter::close() {
	fclose(f_);
}

void DataReader::readInt(int &a) {
	fscanf(f_, "%d\n", &a);
}

void DataPrinter::outputInt(int a) {
	fprintf(f_, "%d\n", a);
}

void DataReader::readSequence(Sequence * &sequence) {
	int length;
	readInt(length);
	if (length == 0) {
		fscanf(f_, "\n");
		sequence = new Sequence("");
	} else {
		char *s = new char[length + 1];
		fscanf(f_, "%s\n", s);
		sequence = new Sequence(s);
	}
}

void DataPrinter::outputSequence(Sequence *sequence) {
	outputInt(sequence->size());
	fprintf(f_, "%s\n", sequence->str().c_str());
}

void DataReader::readEdge(Edge * &edge) {
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
	for (longEdgesMap::iterator it = edges.begin(); it != edges.end(); ++it) {
		if (it->first == it->second->EdgeId) {
			outputInt(it->first);
			outputEdge(it->second);
		}
	}
	Sequence *emptySequence = new Sequence("");
	Edge *emptyEdge = new Edge(emptySequence, emptySequence, 0, 0, 0, 0);
	for (longEdgesMap::iterator it = edges.begin(); it != edges.end(); ++it) {
		if (it->first != it->second->EdgeId) {
			outputInt(it->first);
			emptyEdge->EdgeId = it->second->EdgeId;
			outputEdge(emptyEdge);
		}
	}
	delete emptyEdge;
}

void DataReader::readLongEdgesMap(longEdgesMap &edges) {
	int size;
	readInt(size);
	for (int i = 0; i < size; i++) {
		int id;
		readInt(id);
		Edge *edge;
		readEdge(edge);
		if (id == edge->EdgeId) {
			edges.insert(make_pair(id, edge));
		} else {
			edges.insert(make_pair(id, edges[edge->EdgeId]));
			delete edge;
		}
	}
}

void DataReader::readIntArray(int *array, int length) {
	for (int i = 0; i < length; i++) {
		fscanf(f_, "%d ", array + i);
	}
	fscanf(f_, "\n");
}

void DataPrinter::outputIntArray(int *array, int length) {
	for (int i = 0; i < length; i++) {
		fprintf(f_, "%d ", array[i]);
	}
	fprintf(f_, "\n");
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

void save(char *fileName, PairedGraph &g, longEdgesMap &longEdges,
		int &VertexCount, int EdgeId) {
	DataPrinter dp(fileName);
	dp.outputInt(VertexCount);
	dp.outputInt(EdgeId);
	dp.outputLongEdgesMap(longEdges);
//TODO: FIX!!!
//	dp.outputIntArray(g.inD, MAX_VERT_NUMBER);
//	dp.outputIntArray(g.outD, MAX_VERT_NUMBER);
//	dp.outputIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dp.outputIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dp.close();
}
void load(char *fileName, PairedGraph &g, longEdgesMap &longEdges,
		int &VertexCount, int &EdgeId) {
	DataReader dr(fileName);
	dr.readInt(VertexCount);
	dr.readInt(EdgeId);
	dr.readLongEdgesMap(longEdges);
//TODO: fix;
//	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
//	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
//	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dr.close();
}
