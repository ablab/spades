#include "graphio.hpp"
#include "graphVisualizer.hpp"
#include <stdio.h>

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

void outputLongEdgesThroughGenome(longEdgesMap &longEdges, PairedGraph &graph, int &VertexCount) {
	char Buffer[100];
	char* Genome;
	int GenLength;
	int GenPos;
	int i,t;

	int EdgeNum = 0;

	int bigShift = insertLength+readLength;
	assert(k==l);

	cerr<<"Graph output through genome"<<endl;
	gvis::GraphPrinter<int> g("Paired_ext");

	FILE *infile;
	Genome = (char *) malloc(6000000);
	if ((infile = fopen("data/MG1655-K12_cut.fasta", "r"))
	           == NULL) {
		cerr<<"No such file"<<endl;
		return;
	}
	i=0;
	Genome[i] = fgetc(infile);
	t=0;
	while (t<5) {
	  if (Genome[i]>20) {t=0;i++;}
	  else t++;
	  Genome[i] = fgetc(infile);
	}
	GenLength = i;
	GenPos=0;
	fclose(infile);
	cerr<<"Try to process"<<endl;

	int CurVert = 0;
	while ((graph.inD[CurVert]!=0)||(graph.outD[CurVert]!=1)) CurVert++;
	cerr<<"Start vertex "<<CurVert<<endl;
	while (graph.outD[CurVert]!=0){

		cerr<<"Try to found next edge"<<endl;
		forn(v,graph.outD[CurVert])
		{

			int edgeId = edgeRealId(graph.outputEdges[CurVert][v], longEdges);
			cerr<<"possible edge"<<edgeId<<endl;
			bool goodEdge = true;
			int h=0;
			while(goodEdge&&(h<longEdges[edgeId]->upper->size())){
				//cerr<<" "<<(*(longEdges[edgeId]->upper))[h]<<" =?= "<<Genome[GenPos+h]<<endl;
				if ( nucl((*(longEdges[edgeId]->upper))[h])!=Genome[GenPos+h]) goodEdge=false;
				h++;
			}
			h=0;
			while(goodEdge&&(h<longEdges[edgeId]->lower->size())){
				if (nucl((*(longEdges[edgeId]->lower))[h])!=Genome[GenPos+h+bigShift]) goodEdge=false;
				h++;
			}

			if (goodEdge){
				cerr<<"Edge found"<<endl;
				EdgeNum++;
				sprintf(Buffer, "%i: %i (%i)", EdgeNum, edgeId, longEdges[edgeId]->length);
				g.addEdge(longEdges[edgeId]->FromVertex, longEdges[edgeId]->ToVertex, Buffer);
				cerr << edgeId << " (" << longEdges[edgeId]->length << "):" << endl;
				cerr << longEdges[edgeId]->upper->str() << endl;
				cerr << longEdges[edgeId]->lower->str() << endl;
				CurVert = longEdges[edgeId]->ToVertex;
				GenPos += longEdges[edgeId]->length;
				break;
			}
		}
	}
	g.output();
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

