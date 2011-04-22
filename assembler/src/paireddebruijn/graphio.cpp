#include "graphio.hpp"
#include "graphSimplification.hpp"
#include "graphVisualizer.hpp"
#include <stdio.h>
#include "common.hpp"
#include "pairedGraph.hpp"
#include "iostream"
#include "fstream"

LOGGER("p.graphio");
using namespace paired_assembler;



void codeRead(char *read, char *code) {
	int read_length = strlen(read);
	for (int i = 0; i < read_length; i++) {
		code[i] = codeNucleotide(read[i]);
	}
}

ll extractMer(char *read, int shift, int length) {
	ll res = 0;
	for (int i = 0; i < length; i++) {
		res = res << 2;
		res += codeNucleotide(read[shift + i]);
		if (codeNucleotide( read[shift + i])==-1)
				cerr<<"Extract fault on pos"<<i<<" shift "<<shift;

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

void outputLongEdges(longEdgesMap &longEdges, ostream &os) {
	gvis::GraphPrinter<int> g("Paired_ext", os);
	char Buffer[100];
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second->EdgeId == it->first) {
			sprintf(Buffer, "%i (%i) cov %i", it->first, it->second->length, it->second->coverage);
			//		else sprintf(Buffer,"%i (%i) FAKE now it is %d",it->first, it->second->length,it->second->EdgeId);

			g.addEdge(it->second->FromVertex, it->second->ToVertex, Buffer);
			cerr << it->first << " (" << it->second->length << ") cov "<< it->second->coverage<<":" << endl;
			if (it->second->length < 500) {
				cerr << it->second->upper->str() << endl;
				cerr << it->second->lower->str() << endl;
			}
		}
	}
	g.output();
	cerr<<"Long edge output finished"<<endl;
}

void outputLongEdges(longEdgesMap &longEdges) {
	outputLongEdges(longEdges, cout);
}

void outputLongEdges(longEdgesMap &longEdges, string fileName) {
	ofstream s;
	s.open(fileName.c_str());
	cerr << "Open file " << fileName << endl;
	outputLongEdges(longEdges, s);
	s.close();
}

void outputLongEdges(longEdgesMap &longEdges, PairedGraph &graph, ostream &os) {
	gvis::GraphPrinter<int> g("Paired_ext", os);
	char Buffer[100];
	bool UsedV[200000];
	forn(i,200000)
		UsedV[i] = false;
	pair<int, int> vDist;
	for (longEdgesMap::iterator it = longEdges.begin(); it != longEdges.end(); ++it) {
		if (it->second == NULL) {
			cerr<<"BAD EDGE!!! ID="<<it->first<<endl;
			continue;
		}
		if (it->second->EdgeId == it->first) {
			if (!UsedV[it->second->FromVertex]) {
				vDist = vertexDist(longEdges, graph, it->second->FromVertex);
				sprintf(Buffer, "Vertex_%i (%i, %i)", it->second->FromVertex,
						vDist.first, vDist.second);
				g.addVertex(it->second->FromVertex, Buffer);
				UsedV[it->second->FromVertex] = true;
			}
			if (!UsedV[it->second->ToVertex]) {
				vDist = vertexDist(longEdges, graph, it->second->ToVertex);
				sprintf(Buffer, "Vertex_%i (%i, %i)", it->second->ToVertex,
						vDist.first, vDist.second);
				g.addVertex(it->second->ToVertex, Buffer);
				UsedV[it->second->ToVertex] = true;
			}

			sprintf(Buffer, "%i (%i) cov %i", it->first, it->second->length, it->second->coverage);
			//		else sprintf(Buffer,"%i (%i) FAKE now it is %d",it->first, it->second->length,it->second->EdgeId);

			g.addEdge(it->second->FromVertex, it->second->ToVertex, Buffer);
			cerr << it->first << " (" << it->second->length << ") cov "<< it->second->coverage<<":" << endl;
			if (it->second->length < 500) {
				cerr << it->second->upper->str() << endl;
				cerr << it->second->lower->str() << endl;
			}
		}
	}
	g.output();
}

void outputLongEdges(longEdgesMap &longEdges, PairedGraph &graph) {
	outputLongEdges(longEdges, graph, cout);
}

void outputLongEdges(longEdgesMap &longEdges, PairedGraph &graph,
		string fileName) {
	ofstream s;
	s.open(fileName.c_str());
	outputLongEdges(longEdges, graph, s);
	s.close();
}

Sequence readGenome(istream &is) {
	SequenceBuilder sb;
	string buffer;
	while(!is.eof()){
		is >> buffer;
		sb.append(Sequence(buffer));
	}
	return sb.BuildSequence();
}

Sequence readGenomeFromFile(const string &fileName) {
	ifstream is;
	is.open(fileName.c_str());
	Sequence result(readGenome(is));
	is.close();
	return result;
}

int findStartVertex(PairedGraph &graph, Sequence &genome, int position = 0) {
	int result = -1;
//	cerr<<"findStartVertex"<<endl;
	for (int i = 0; i < graph.VertexCount; i++) {
//		if (graph.degrees[i][0] == 0 && graph.degrees[i][1] == 1) {
		for (int edgeNum = 0; edgeNum<graph.degrees[i][1]; edgeNum++){
	//			cerr<<"SEQ VS GEN"<<endl;
			Sequence* tmp_seq = graph.longEdges[graph.edgeIds[i][edgeNum][OUT_EDGE]]->upper;
//			cerr<<"Seq "<<tmp_seq->str()<<endl;
//			cerr<<"Gen "<<(genome.Subseq(0,tmp_seq->size())).str()<<endl;
			if (genome.Subseq(position,tmp_seq->size()+position)== *tmp_seq){
				if (result >= 0) {
					cerr << "Ambigious start point for threading!" << endl;
					return result;
				}
				result = i;
			}
		}
	}
	return result;
}

bool checkEdge(Edge *nextEdge, int genPos, Sequence &genome) {
	for (size_t i = 0; i < nextEdge->upper->size(); i++)
		if (nextEdge->upper->operator [](i) != genome[genPos + i]
//				|| nextEdge->lower->operator [](i) != genome[genPos + i
//						+ readLength + insertLength]
						)
			return false;
	return true;
}

string createEdgeLabel(int edgeNum, int edgeId, int length) {
	stringstream ss;
	ss << edgeNum << ":" << edgeId << " (" << length << ")";
	return ss.str();
}

int moveThroughEdge(gvis::GraphPrinter<int> &g, PairedGraph &graph,
		Edge *nextEdge, int edgeNum, int genPos) {
	cerr << "Edge found" << endl;
	edgeNum++;
	string label(createEdgeLabel(edgeNum, nextEdge->EdgeId, nextEdge->length));
	g.addEdge(nextEdge->FromVertex, nextEdge->ToVertex, label, "red");
/*	cerr << nextEdge->EdgeId << " (" << nextEdge->length << "):" << endl;
	if (graph.longEdges[nextEdge->EdgeId]->length < 500) {
		cerr << nextEdge->upper->str() << endl;
		cerr << nextEdge->lower->str() << endl;
	}
*/	genPos += nextEdge->length;
	return nextEdge->ToVertex;
}

Edge *chooseNextEdge(int currentVertex, int genPos, Sequence &genome,
		PairedGraph &graph) {
	for (int v = 0; v < graph.degrees[currentVertex][1]; v++) {
		int notRealId = graph.edgeIds[currentVertex][v][OUT_EDGE];
		int edgeId = edgeRealId(notRealId, graph.longEdges);
		Edge *nextEdge = graph.longEdges[edgeId];
		cerr << "possible edge" << edgeId << endl;
		if (checkEdge(nextEdge, genPos, genome)) {
			return nextEdge;
		}
	}
	return NULL;
}
#define MAX_GAP 2000
void outputLongEdgesThroughGenome(PairedGraph &graph, ostream &os) {
	assert(k==l);
	cerr << "Graph output through genome" << endl;
	gvis::GraphPrinter<int> g("Paired_ext", os);
	Sequence genome(readGenomeFromFile("data/input/MG1655-K12_cut.fasta"));
	cerr << "Try to process" << endl;
	int gap = 0;
	int edgeNum = 0;
	int genPos = 0;
	int currentVertex = findStartVertex(graph, genome);
	cerr << "Start vertex " << currentVertex << endl;
//	while (graph.degrees[currentVertex][1] != 0) {
	while (1) {
		cerr << "Try to found next edge" << endl;
		Edge *nextEdge = chooseNextEdge(currentVertex, genPos, genome, graph);
		if (nextEdge != NULL) {
			currentVertex
					= moveThroughEdge(g, graph, nextEdge, edgeNum, genPos);
			edgeNum++;
			genPos += nextEdge->length;
		} else {
			gap = 0;
			cerr<<"Gap sequence from pos "<<genPos<<":"<<endl;
			while (((currentVertex = findStartVertex(graph, genome, gap+genPos))==-1)&&(gap<MAX_GAP)){
				cerr<<nucl(genome[gap+genPos]);
				gap++;
			}
			cerr<<endl;
			if (gap >= MAX_GAP) {
				cerr << "BAD GRAPH. I can not cover all genome." << endl;
				break;
			}
			else {
				genPos+=gap;
			}
		}
	}
	cerr<<"Go trough the graph finished on position "<<genPos+k-1<<endl;
	g.output();
}

void outputLongEdgesThroughGenome(PairedGraph &graph) {
	outputLongEdgesThroughGenome(graph, cout);
}

void outputLongEdgesThroughGenome(PairedGraph &graph, string fileName) {
	ofstream s;
	s.open(fileName.c_str());
	outputLongEdgesThroughGenome(graph, s);
	s.close();
}

DataReader::DataReader(char *fileName) {
	f_ = fopen(fileName, "r");
	DEBUG("DataReader " << fileName <<" created");
	assert(f_ != NULL);
}

DataReader::DataReader(const char *fileName) {
	f_ = fopen(fileName, "r");
	DEBUG("DataReader " << fileName <<" created");
	assert(f_ != NULL);
}

DataPrinter::DataPrinter(char *fileName) {
	f_ = fopen(fileName, "w");
	assert(f_ != NULL);
}

DataPrinter::DataPrinter(const char *fileName) {
	f_ = fopen(fileName, "w");
	assert(f_ != NULL);
}

void DataReader::close() {
	fclose(f_);
}

void DataPrinter::close() {
	fclose(f_);
}

int DataReader::read(int &a) {
	return fscanf(f_, "%d\n", &a);
}

void DataPrinter::output(int a) {
	fprintf(f_, "%d\n", a);
}

int DataReader::read(long long &a) {
	return fscanf(f_, "%lld\n", &a);
}

void DataPrinter::output(long long a) {
	fprintf(f_, "%lld\n", a);
}

int DataReader::read(Sequence * &sequence) {
	int length;
	if (!read(length)) {
		return 0;
	} else {
		if (length == 0) {
			fscanf(f_, "\n");
			sequence = new Sequence("");
		} else {
			char *s = new char[length + 1];
			fscanf(f_, "%s\n", s);
			sequence = new Sequence(s);
		}
		return 1;
	}
}

void DataPrinter::output(Sequence *sequence) {
	output((int) sequence->size());
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

int DataReader::read(Edge * &edge) {
	int from, to, len, id, cov;
	Sequence *up, *low;
	int read_res = 0;
	read_res += read(id);
	if (id == -1) return 0;
	DEBUG(id);
	read_res += read(from);
	read_res += read(to);
	read_res += read(len);
	read_res += read(cov);
	read_res += read(up);
	read_res += read(low);
	edge = new Edge(up, low, from, to, len, id, cov);
	DEBUG("edge"<<id << "loaded");
	return (read_res/7);
}

void DataPrinter::output(Edge *edge) {
	output(edge->EdgeId);
	output(edge->FromVertex);
	output(edge->ToVertex);
	output(edge->length);
	output(edge->coverage);
	output(edge->upper);
	output(edge->lower);
}

void DataPrinter::outputLongEdgesMap(longEdgesMap &edges) {
	INFO("Saving long edges");
	DEBUG(edges.size());
	int size = edges.size();
	for (longEdgesMap::iterator it = edges.begin(); it != edges.end(); ++it) {
		if (it->first == it->second->EdgeId) {
			output(it->first);
			output(it->second);
			DEBUG("Edge outputed" );
		}
	}
	DEBUG("Normal edges outputed");
	Sequence *emptySequence = new Sequence("");
	Edge *emptyEdge = new Edge(emptySequence, emptySequence, 0, 0, 0, 0);
	for (longEdgesMap::iterator it = edges.begin(); it != edges.end(); ++it) {
		if (it->first != it->second->EdgeId) {
			output(it->first);
			emptyEdge->EdgeId = it->second->EdgeId;
			output(emptyEdge);
		}
	}
	DEBUG("fakeEdges outputed");
	delete emptyEdge;
	output(-1);
	output(size);
}

void DataReader::readLongEdgesMap(longEdgesMap &edges) {
	int size = 0;
	int id;
	while (1) {
		assert(read(id));
		Edge *edge;
		DEBUG(id);
		if( id == -1 || !(read(edge)))
			break;
		else {
			size++;
			if (id == edge->EdgeId) {
				edges.insert(make_pair(id, edge));
			} else {
				edges.insert(make_pair(id, edges[edge->EdgeId]));
				delete edge;
			}
		}
	}
	read(id);
	DEBUG(id);
	assert(size == id);
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

void save(DataPrinter dp, Edge *e){
	dp.output(e->EdgeId);
	dp.output(e);
//	dp.close();
}
void save(char *fileName, PairedGraph &g, longEdgesMap &longEdges,
		int &VertexCount, int EdgeId) {
	DataPrinter dp(fileName);
	dp.outputLongEdgesMap(longEdges);
	dp.output(VertexCount);
	dp.output(EdgeId);
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
	dr.readLongEdgesMap(longEdges);
	dr.read(VertexCount);
	dr.read(EdgeId);
	//TODO: fix;
	//	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
	//	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
	//	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	//	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dr.close();
}

void save(char *fileName, PairedGraph &g) {
	INFO("Saving graph");
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
void save(string fileName, PairedGraph &g) {
	INFO("Saving graph");
	DataPrinter dp(fileName.c_str());
	dp.outputLongEdgesMap(g.longEdges);
	dp.output(g.VertexCount);
	dp.output(g.EdgeId);
	//TODO: FIX!!!
	//	dp.outputIntArray(g.inD, MAX_VERT_NUMBER);
	//	dp.outputIntArray(g.outD, MAX_VERT_NUMBER);
	//	dp.outputIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	//	dp.outputIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dp.close();
}
void load(char *fileName, PairedGraph &g) {
	INFO("Loading graph");
	DataReader dr(fileName);
	dr.readLongEdgesMap(g.longEdges);
	dr.read(g.VertexCount);
	dr.read(g.EdgeId);
	//TODO: fix;
	//	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
	//	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
	//	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	//	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dr.close();
}
void load(string fileName, PairedGraph &g) {
	INFO("Loading graph");
	DataReader dr(fileName.c_str());
	dr.readLongEdgesMap(g.longEdges);
	dr.read(g.VertexCount);
	dr.read(g.EdgeId);
	//TODO: fix;
	//	dr.readIntArray(g.inD, MAX_VERT_NUMBER);
	//	dr.readIntArray(g.outD, MAX_VERT_NUMBER);
	//	dr.readIntArray((int*) g.outputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	//	dr.readIntArray((int*) g.inputEdges, MAX_VERT_NUMBER, MAX_DEGREE);
	dr.close();
}

void outputVertexKmers(edgesMap &edges){

	FILE *fkmers = fopen((folder+"kmers.txt").c_str(), "w");
	for (edgesMap::iterator iter = edges.begin(); iter != edges.end();++iter) {
		ll kmer = iter->fi;
		fprintf(fkmers,"%lld %s\n", kmer, decompress(kmer, k));
	}
	fclose(fkmers);

}
