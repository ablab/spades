#include "graphio.hpp"
#include "graphSimplification.hpp"
#include "graphVisualizer.hpp"
#include <stdio.h>
#include "common.hpp"
#include "pairedGraph.hpp"
#include "iostream"
#include "fstream"

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

void outputLongEdges(longEdgesMap &longEdges, string fileName) {
	gvis::GraphPrinter<int> *g;
	ofstream s;
	if(fileName == "") {
		g = new gvis::GraphPrinter<int>("Paired_ext", cout);
	} else {
		s.open(fileName.c_str());
		cerr<<"Open file "<<fileName<<endl;
		g = new gvis::GraphPrinter<int>("Paired_ext", s) ;
	}
	char Buffer[100];
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			sprintf(Buffer, "%i (%i)", it->first, it->second->length);
			//		else sprintf(Buffer,"%i (%i) FAKE now it is %d",it->first, it->second->length,it->second->EdgeId);

			g->addEdge(it->second->FromVertex, it->second->ToVertex, Buffer);
			cerr << it->first << " (" << it->second->length << "):" << endl;
			if (it->second->length < 500) {
				cerr << it->second->upper->str() << endl;
				cerr << it->second->lower->str() << endl;
			}
		}
	}
	g->output();
	if(fileName != "") {
		s.close();
	}
	delete g;
}

void outputLongEdges(longEdgesMap &longEdges, PairedGraph &graph, string fileName) {
	gvis::GraphPrinter<int> *g;
	ofstream s;
	if(fileName == "") {
		g = new gvis::GraphPrinter<int>("Paired_ext", cout);
	} else {
		s.open(fileName.c_str());
		g = new gvis::GraphPrinter<int>("Paired_ext", s) ;
	}
	char Buffer[100];
	bool UsedV[20000];
	forn(i,20000) UsedV[i] = false;
	pair<int,int> vDist;
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			if (!UsedV[it->second->FromVertex]){
				vDist = vertexDist(longEdges, graph,it->second->FromVertex);
				sprintf(Buffer, "Vertex_%i (%i, %i)", it->second->FromVertex, vDist.first, vDist.second );
				g->addVertex(it->second->FromVertex, Buffer);
				UsedV[it->second->FromVertex]=true;
			}
			if (!UsedV[it->second->ToVertex]){
				vDist = vertexDist(longEdges, graph,it->second->ToVertex);
				sprintf(Buffer, "Vertex_%i (%i, %i)", it->second->ToVertex, vDist.first, vDist.second );
				g->addVertex(it->second->ToVertex, Buffer);
				UsedV[it->second->ToVertex]=true;
			}

			sprintf(Buffer, "%i (%i)", it->first, it->second->length);
			//		else sprintf(Buffer,"%i (%i) FAKE now it is %d",it->first, it->second->length,it->second->EdgeId);
			g->addEdge(it->second->FromVertex, it->second->ToVertex, Buffer);
			cerr << it->first << " (" << it->second->length << "):" << endl;
			if (it->second->length < 500)
			{
			cerr << it->second->upper->str() << endl;
			cerr << it->second->lower->str() << endl;
			}
		}
	}
	g->output();
	if(fileName != "") {
		s.close();
	}
	delete g;

}


void outputLongEdgesThroughGenome(PairedGraph &graph, string fileName) {
	gvis::GraphPrinter<int> *g;
	ofstream s;
	if(fileName == "") {
		g = new gvis::GraphPrinter<int>("Paired_ext", cout);
	} else {
		s.open(fileName.c_str());
		g = new gvis::GraphPrinter<int>("Paired_ext", s) ;
	}
	char Buffer[100];
	char* Genome;
	int GenLength;
	int GenPos;
	int i, t;

	int EdgeNum = 0;

	int bigShift = insertLength + readLength;
	assert(k==l);

	cerr << "Graph output through genome" << endl;

	FILE *infile;
	Genome = (char *) malloc(6000000);
	if ((infile = fopen("data/MG1655-K12_cut.fasta", "r")) == NULL) {
		cerr << "Can not find genome file data/MG1655-K12_cut.fasta" << endl;
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
			int edgeId = edgeRealId(graph.edgeIds[CurVert][v][OUT_EDGE], graph.longEdges);
			cerr << "possible edge" << edgeId << endl;
			bool goodEdge = true;
			int h = 0;
			while (goodEdge && (h < graph.longEdges[edgeId]->upper->size())) {
				//cerr<<" "<<(*(graph.longEdges[edgeId]->upper))[h]<<" =?= "<<Genome[GenPos+h]<<endl;
				if (nucl((*(graph.longEdges[edgeId]->upper))[h])
						!= Genome[GenPos + h])
					goodEdge = false;
				h++;
			}
			h = 0;
			while (goodEdge && (h < graph.longEdges[edgeId]->lower->size())) {
				if (nucl((*(graph.longEdges[edgeId]->lower))[h]) != Genome[GenPos + h
						+ bigShift])
					goodEdge = false;
				h++;
			}

			if (goodEdge) {
				cerr << "Edge found" << endl;
				EdgeNum++;
				sprintf(Buffer, "%i: %i (%i)", EdgeNum, edgeId,
						graph.longEdges[edgeId]->length);
				g->addEdge(graph.longEdges[edgeId]->FromVertex,
						graph.longEdges[edgeId]->ToVertex, Buffer);
				cerr << edgeId << " (" << graph.longEdges[edgeId]->length << "):"
						<< endl;
				if (graph.longEdges[edgeId]->length < 500) {
					cerr << graph.longEdges[edgeId]->upper->str() << endl;
					cerr << graph.longEdges[edgeId]->lower->str() << endl;
				}
				CurVert = graph.longEdges[edgeId]->ToVertex;
				GenPos += graph.longEdges[edgeId]->length;
				NoEdge = false;
				break;
			}
		}
		if (NoEdge) {
			cerr<<"BAD GRAPH. I can not cover all genome"<<endl;
			break;
		}
	}
	g->output();
	if(fileName != "") {
		s.close();
	}
	delete g;

}

DataReader::DataReader(char *fileName) {
	f_ = fopen(fileName, "r");
	assert(f_ != NULL);
}

DataPrinter::DataPrinter(char *fileName) {
	f_ = fopen(fileName, "w");
	assert(f_ != NULL);
}

void DataReader::close() {
	fclose(f_);
}

void DataPrinter::close() {
	fclose(f_);
}

void DataReader::read(int &a) {
	assert(fscanf(f_, "%d\n", &a)==1);
}

void DataPrinter::output(int a) {
	fprintf(f_, "%d\n", a);
}

void DataReader::read(long long &a) {
	fscanf(f_, "%lld\n", &a);
}

void DataPrinter::output(long long a) {
	fprintf(f_, "%lld\n", a);
}

void DataReader::read(Sequence * &sequence) {
	int length;
	read(length);
	if (length == 0) {
		fscanf(f_, "\n");
		sequence = new Sequence("");
	} else {
		char *s = new char[length + 1];
		fscanf(f_, "%s\n", s);
		sequence = new Sequence(s);
	}
}

void DataPrinter::output(Sequence *sequence) {
	output((int)sequence->size());
	fprintf(f_, "%s\n", sequence->str().c_str());
}

void DataReader::read(VertexPrototype * &v) {
	int id;
	Sequence *lower;
	bool b;
	read(id);
	read(lower);
	int tmpInt;
	read(tmpInt);
	b = tmpInt;
	v = new VertexPrototype(lower, id);
	v->used = b;
}

void DataPrinter::output(VertexPrototype *v) {
	output(v->VertexId);
	output(v->lower);
	output(v->used);
}

void DataReader::read(Edge * &edge) {
	int from, to, len, id;
	Sequence *up, *low;
	read(id);
	read(from);
	read(to);
	read(len);
	read(up);
	read(low);
	edge = new Edge(up, low, from, to, len, id);
}

void DataPrinter::output(Edge *edge) {
	output(edge->EdgeId);
	output(edge->FromVertex);
	output(edge->ToVertex);
	output(edge->length);
	output(edge->upper);
	output(edge->lower);
}

void DataPrinter::outputLongEdgesMap(longEdgesMap &edges) {
	output((int)edges.size());
	for (longEdgesMap::iterator it = edges.begin(); it != edges.end(); ++it) {
		if (it->first == it->second->EdgeId) {
			output(it->first);
			output(it->second);
		}
	}
	Sequence *emptySequence = new Sequence("");
	Edge *emptyEdge = new Edge(emptySequence, emptySequence, 0, 0, 0, 0);
	for (longEdgesMap::iterator it = edges.begin(); it != edges.end(); ++it) {
		if (it->first != it->second->EdgeId) {
			output(it->first);
			emptyEdge->EdgeId = it->second->EdgeId;
			output(emptyEdge);
		}
	}
	delete emptyEdge;
}

void DataReader::readLongEdgesMap(longEdgesMap &edges) {
	int size;
	read(size);
	for (int i = 0; i < size; i++) {
		int id;
		read(id);
		Edge *edge;
		read(edge);
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
	dp.output(VertexCount);
	dp.output(EdgeId);
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
	dr.read(VertexCount);
	dr.read(EdgeId);
	dr.readLongEdgesMap(longEdges);
//TODO: fix;
//	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
//	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
//	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dr.close();
}

void save(char *fileName, PairedGraph &g) {
	DataPrinter dp(fileName);
	dp.output(g.VertexCount);
	dp.output(g.EdgeId);
	dp.outputLongEdgesMap(g.longEdges);
//TODO: FIX!!!
//	dp.outputIntArray(g.inD, MAX_VERT_NUMBER);
//	dp.outputIntArray(g.outD, MAX_VERT_NUMBER);
//	dp.outputIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dp.outputIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dp.close();
}
void load(char *fileName, PairedGraph &g) {
	DataReader dr(fileName);
	dr.read(g.VertexCount);
	dr.read(g.EdgeId);
	dr.readLongEdgesMap(g.longEdges);
//TODO: fix;
//	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
//	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
//	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
//	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dr.close();
}
