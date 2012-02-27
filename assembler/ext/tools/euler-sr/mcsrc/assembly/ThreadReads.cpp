/***************************************************************************
 * Title:          ThreadReads.cpp 
 * Author:         Dumitru Brinza dima@cs.ucsd.edu
 * Created:        2008
 * Last modified:  11/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "DeBruijnGraph.h"
#include "IntervalGraph.h"

#include "DNASequence.h"
#include "SeqReader.h"
#include "utils.h"
#include "string.h"
#include "IntegralTupleStatic.h"


using namespace std;

//---------------------------------------------------------------------------------
void PrintUsage() {
cout << "usage: threadReads graphBase readsPrefix originalReads outputFile [vertexSize] [hammingDistance] [outputEdges]" << endl << endl;

cout << "graphBase       -- name of the files containing graph information                " << endl;
cout << "readsPrefix     -- file containing fixed by voting prefixes of reads             " << endl;
cout << "originalReads   -- file containing original set of reads. Read ID should be equal" << endl;
cout << "                   to index in the titles of reads from readsPrefix              " << endl;
cout << "outputFile      -- name of the file with outputted threaded reads                " << endl;
cout << "vertexSize      -- size of the vertex in the graph, default 30                   " << endl;
cout << "hammingDistance -- use Hamming distance when choosing the clothest to read path, " << endl;
cout << "                   value 0 - edit distance, 1 - Hamming, default is 0            " << endl;
cout << "outputEdges     -- append graph edges to the outputFile, default is 1, 0 means no" << endl;

}

#define MAXSTRING 2000		// Maximum size for reading in a string

//---------------------------------------------------------------------------------
ssize_t  PrintEdges(IntervalGraph &graph,std::ofstream &contigOut);
void TraversEdges(IntervalGraph &graph, DNASequenceList &reads, vector<string>& threadedReads,ssize_t first_ndex, ssize_t last_index);
ssize_t  StorePathSequences(IntervalGraph &graph, ssize_t edgeIndex,ssize_t edgePos, std::string curSequence,
                        ssize_t searchLength, std::vector<std::string> &sequences, string trace_path, 
                        vector<ssize_t> segment, vector<ssize_t> &resulted_segment,ssize_t level=0);
void Thread(IntervalGraph &graph,ssize_t edgeIndex, ssize_t intvIndex, DNASequence &read, vector<string>& threadedReads,
                        bool HammingDistance=false);
void FindBestMatch(DNASequence &read, std::vector<std::string> &pathSequences, ssize_t & id_best, 
                        ssize_t & id_best2, ssize_t & delta_best, ssize_t &delta_best2, bool HammingDistance);
ssize_t  EditDistance(char A[], char B[], ssize_t match_penalty=1, ssize_t insert_penalty=1, ssize_t delete_penalty=1);
void ExtendReads(DNASequenceList &prefix_reads, vector<string> &original_reads);
bool ReadOneSeq(ifstream & infile, string & title, string & sequence);
ssize_t  ReadFast(string prefixFile, string originalFile, string outputFile);
//---------------------------------------------------------------------------------

ssize_t  VERTEX_SIZE  = 30;
bool HAMMING_DIST = false;
bool OUTPUT_EDGES = 1;

SimpleSequence seq;
DNASequence read2;
//---------------------------------------------------------------------------------
int main(int argc, char* argv[]) {

  if (argc < 5) {
                 PrintUsage();
                 exit(1);
                }
  
  string graphBase = argv[1];
  string prefixFile = argv[2];
  string originalFile = argv[3];
  string outputFile = argv[4];
  if(argc>5)VERTEX_SIZE = atoi(argv[5]);
  if(argc>6)HAMMING_DIST = atoi(argv[6]);
  if(argc>7)OUTPUT_EDGES = atoi(argv[7]);

  cout << "VERTEX SIZE = " << VERTEX_SIZE << endl;
  cout << "HAMMING DIST = " << HAMMING_DIST << endl;
  cout << "OUTPUT EDGES = " << OUTPUT_EDGES << endl;

  string graphFile = graphBase + ".bgraph";
  string intvFile  = graphBase + ".intv";
  string pathFile  = graphBase + ".path";
  string edgeFile  = graphBase + ".edge";

  IntervalGraph graph;
  graph.vertexSize = VERTEX_SIZE;
  read2.length = 0;
  read2.seq = new unsigned char[MAXSTRING];
	
  cout << "Read graph ..." << endl;	
  graph.ReadIntervalGraph(graphFile, intvFile, pathFile);
  cout << "Read edges ..." << endl;
  ReadSequences(edgeFile, graph.edges);

  cout << "Read file with prefixes and original reads, output concatenation ... " << endl;	
  ssize_t nReads = ReadFast(prefixFile, originalFile, outputFile + ".tmp" );
  cout << "Total # reads = " << nReads << endl;

  cout << "Thread and output reads ... " << endl;	
  
  ifstream concatfile;
  openck(outputFile + ".tmp", concatfile, std::ios::in);

  std::ofstream contigOut;
  openck(outputFile, contigOut, std::ios::out);

  string concat_sequence = "";
  string concat_title="",concat_title2="";
    
  //UNUSED+// ssize_t count=0;
  ssize_t first_index=0,last_index=-1 ;
  DNASequenceList concat_reads;
  DNASequence sequence;
  
  vector<string> threadedReads;

  ReadOneSeq(concatfile, concat_title2, concat_sequence);
  while(ReadOneSeq(concatfile, concat_title, concat_sequence))
  {
    if(++last_index%100000==0&&last_index!=0){

     vector<string> threadedReads;
     TraversEdges(graph,concat_reads,threadedReads,first_index, last_index);
     cout << last_index*100/nReads << "% " << threadedReads.size() << endl;	
     for(ssize_t i=0;i< threadedReads.size();i++) contigOut << threadedReads[i] << endl;
     first_index=last_index;
     concat_reads.clear();
    
    }
    
   sequence.namestr = concat_title2.substr(1);
   sequence.Reset(concat_sequence.length());
   memcpy(&sequence.seq[0], concat_sequence.c_str(), concat_sequence.length());
   sequence.length = concat_sequence.length();
   concat_reads.push_back(sequence);
   
   concat_title2=concat_title;
   }
  
  concatfile.close();
  
  if(OUTPUT_EDGES)
  {
  cout << "Output edges ... " << endl;	
  PrintEdges(graph,contigOut);
  }

 contigOut.close();
 return 0;
}

//---------------------------------------------------------------------------------
void ExtendReads(string &prefix_read, string &original_read)
 {
   for(ssize_t j=0;j<prefix_read.length();j++) prefix_read[j] = toupper(prefix_read[j]);
   for(ssize_t j=0;j<original_read.length();j++) original_read[j] = toupper(original_read[j]);
         
   ssize_t pos = original_read.find(prefix_read.substr(0,5));
   if(pos<0) pos = original_read.find(prefix_read.substr(5,5),5)-5;
   if(pos<0) pos = original_read.find(prefix_read.substr(10,5),10)-10;
           
   if(pos<0)
   {
    cout << " errors possible " << endl;
    pos=0;
   }
           
   ssize_t extension = original_read.length() - prefix_read.length() - pos;
   if(extension>0) prefix_read +=  original_read.substr(prefix_read.length()+pos,extension);
 }
//---------------------------------------------------------------------------------
void TraversEdges(IntervalGraph &graph, DNASequenceList &reads, vector<string>& threadedReads,ssize_t first_index, ssize_t last_index)
{
	ssize_t edgeIndex;
	ssize_t readIndex;
	char *edgeSequence;
	ssize_t edgeSequencLength;
	ssize_t intv;
	ssize_t step = 0;
	
	for (edgeIndex = 0 ; edgeIndex < graph.edges.size(); edgeIndex++ ){
	         if(step++>graph.edges.size()/50){cout<<".";cout.flush();step=0;}
		// quick access to the edge
		edgeSequence = (char *) graph.edges[edgeIndex].seq.seq;
		edgeSequencLength = graph.edges[edgeIndex].seq.length;

		// Now, to try and thread reads through the graph,
		// find out what reads are mapped to this edge.
		// These are stored in the interval list.

		for (intv = 0; intv < (*graph.edges[edgeIndex].intervals).size(); intv++) {
			if ((*graph.edges[edgeIndex].intervals)[intv].readPos == 0) {
				readIndex = (*graph.edges[edgeIndex].intervals)[intv].read;
				if(readIndex % 2 == 0 && readIndex/2>=first_index && readIndex/2<last_index) // Skip complement sequence
				Thread(graph, edgeIndex, intv, reads[readIndex/2-first_index], threadedReads,HAMMING_DIST);
			}
		}

	}
	cout << endl;

}
//-------------------------------------------------------------------------------------
string itos(ssize_t i)	// convert int to string
{
		stringstream s;
		s << i;
		return s.str();
}
//-------------------------------------------------------------------------------------
void Thread(IntervalGraph &graph,ssize_t edgeIndex, ssize_t intvIndex, DNASequence &read, vector<string>& threadedReads, bool HammingDistance) {

	std::vector<std::string> pathSequences;
	
	//Store all possible path for the read
	vector<ssize_t>   segment,resulted_segment;
	StorePathSequences(graph, edgeIndex, (*graph.edges[edgeIndex].intervals)[intvIndex].edgePos, "", read.length, pathSequences,"",segment,resulted_segment,0);         
          
          //Find the best (closest) path
	//UNUSED+// ssize_t num;
	ssize_t id_best, id_best2, delta_best, delta_best2,  seqlength;
	FindBestMatch(read, pathSequences,id_best,id_best2,delta_best,delta_best2,HammingDistance);
	
	// Output reads which have at least 2 alternative paths
	if(id_best>-1) 
	{
           if(!HammingDistance)
	  {
	    if(id_best2>-1 && delta_best<read.length*0.2)
	    {
	    threadedReads.push_back(">"+read.namestr);
	    threadedReads.push_back(pathSequences[id_best]);
	    }
  	  }
	else  
	  { 
               seqlength = pathSequences[id_best].length();
               if(seqlength>=read.length)seqlength = read.length;

               if(delta_best>0)
               {
               ssize_t total=0, errors=0, catastrophy_begin=seqlength;
               for(ssize_t i=seqlength-1;i>VERTEX_SIZE;i--)
               {
                total++;
                if(pathSequences[id_best][i]!=read.seq[i]){errors++;}
                if((double)errors/total>=0.50){catastrophy_begin=i-1;}
               }
               
               seqlength = catastrophy_begin;
               read.length = seqlength;
               
               if(id_best2>-1&&delta_best2-delta_best<=3)
               {
               string bestSeq = pathSequences[id_best].substr(0,seqlength);
               StorePathSequences(graph, edgeIndex, (*graph.edges[edgeIndex].intervals)[intvIndex].edgePos, "", seqlength, pathSequences,pathSequences[id_best],segment,resulted_segment,1);
               
               read2.length = seqlength;
               for(ssize_t i=0;i<read2.length;i++) read2.seq[i] = bestSeq[i];
               FindBestMatch(read2, pathSequences,id_best,id_best2,delta_best,delta_best2,HammingDistance);
               if(id_best2>-1&&delta_best2-delta_best<=3)
               {
                bool was_cut = false;
                if(resulted_segment.size()>0)
                {
                 for(ssize_t i=seqlength-(resulted_segment[resulted_segment.size()-1]-VERTEX_SIZE);i<seqlength&&pathSequences[id_best2].length();i++)
                 {         
	        if(pathSequences[id_best][i]!=pathSequences[id_best2][i]){
	        was_cut=true;
	        seqlength -=resulted_segment[resulted_segment.size()-1]-VERTEX_SIZE;
	        resulted_segment.erase(resulted_segment.end()-1);
	        break;
	        }
	       }
	      }
               if(!was_cut)  seqlength = VERTEX_SIZE;
               }}
               }
 
                threadedReads.push_back(">"+read.namestr);
	      threadedReads.push_back(pathSequences[id_best].substr(0,seqlength));
          }}
}
//-------------------------------------------------------------------------------------------
ssize_t StorePathSequences(IntervalGraph &graph, ssize_t edgeIndex,ssize_t edgePos, std::string curSequence,ssize_t searchLength, std::vector<std::string> &sequences, string trace_path, 
											 vector<ssize_t> segment, vector<ssize_t> &resulted_segment,ssize_t level) {
	std::string newSequence;
	if (searchLength == 0) 
          return sequences.size();
	//Check if read ends on the current edge  
	if (edgePos + searchLength <= graph.edges[edgeIndex].length ) {
		newSequence = "";
		for(ssize_t i=0;i<searchLength;i++) {newSequence +=" "; newSequence[i] = graph.edges[edgeIndex].seq.seq[edgePos+i];}

                    if(level!=0)
                    {
                     curSequence += newSequence;
                     ssize_t l1 = trace_path.length();
                     ssize_t l2 = curSequence.length();
                     ssize_t lmin = (l1>l2?l2:l1);
                     if(curSequence.substr(0,lmin)==trace_path.substr(0,lmin))
                     {
                     for(ssize_t i=0;i<segment.size();i++)  resulted_segment.push_back(segment[i]);
                     resulted_segment.push_back(newSequence.length());
                     }
                     }
                     sequences.push_back(curSequence + newSequence);


	}
	else {
		// Not done searching graph for sequence of length 'searchLength'
		newSequence = "";

		ssize_t dest;
		dest = graph.edges[edgeIndex].dest;
		ssize_t destVertexSize = VERTEX_SIZE;
		for(ssize_t i=0;i<graph.edges[edgeIndex].length - destVertexSize - edgePos;i++) {newSequence +=" "; newSequence[i] = graph.edges[edgeIndex].seq.seq[edgePos+i];}

		// Append the substring
		curSequence = curSequence + newSequence;
			
		if(level!=0)
			{
				segment.push_back(newSequence.length());
			}
                                         
		searchLength -= (graph.edges[edgeIndex].length - destVertexSize - edgePos);
		
		ssize_t outEdgeIndex, outEdge, branchNumber = 0;
		for (outEdgeIndex = graph.vertices[dest].FirstOut();
				 outEdgeIndex < graph.vertices[dest].EndOut();
				 outEdgeIndex = graph.vertices[dest].NextOut(outEdgeIndex)) {
				 outEdge = graph.vertices[dest].out[outEdgeIndex];
			
			branchNumber++;
			StorePathSequences(graph, outEdge, 0, curSequence, searchLength, sequences,trace_path,segment,resulted_segment,(level>0?level+1:0));
			}
		if(branchNumber==0&&curSequence.length()>=VERTEX_SIZE) sequences.push_back(curSequence);	
		
	}
	return sequences.size();
}

//----------------------------------------------------------------------
void FindBestMatch(DNASequence &read, std::vector<std::string> &pathSequences, ssize_t & id_best, ssize_t & id_best2, ssize_t & delta_best, ssize_t &delta_best2, bool HammingDistance)
{
 id_best = -1;
 id_best2 = -1;
 delta_best = MAXSTRING;
 delta_best2 = MAXSTRING;

 
 if(pathSequences.size()==1)
  {
   id_best = 0;delta_best = 0;
   return;
  }

 char b[read.length];
 memcpy(b,read.seq,read.length);
 
	for(ssize_t i=0;i<pathSequences.size();i++)
	{
	 ssize_t n_miss=0;
	 if(HammingDistance){
	  //Hamming Distance
	    for(ssize_t x=0;(x<pathSequences[i].length()&&x<100/*&&x<read.length*/);x++) if(read.seq[x]!=pathSequences[i][x]) n_miss++;
	  }
	  else
	  {
	  //Edit distance 
	    ssize_t lon = pathSequences[i].length();
	    if(lon>read.length)lon = read.length;
	    char a[lon];
	    memcpy(a,pathSequences[i].c_str(),lon);
              n_miss = EditDistance(a,b,1,2,2);// penalty for mismatch, insertion, deletion
            }
	 if(n_miss<=delta_best) {delta_best2 = delta_best; delta_best = n_miss; id_best2 = id_best; id_best = i;}
	 else
 	 if(n_miss<=delta_best2) {delta_best2 = n_miss; id_best2 = i;}
          }
}

//-----------------------------------------------------------------------

ssize_t PrintEdges(IntervalGraph &graph,std::ofstream &contigOut) {

	std::stringstream titleStrm;
	for (ssize_t e = 0; e < graph.edges.size(); e++ ) {
		{
			titleStrm.str(""); 
			titleStrm << graph.edges[e].src << " -> " << graph.edges[e].dest << " (" << e << ")";
			graph.edges[e].seq.PrintSeq(contigOut, titleStrm.str());
		}
	}
	return 0;
}
//--------------------------------------------------------------------------------------
bool ReadOneSeq(ifstream & infile, string & title, string & sequence)
{
  string line;
  sequence = "";
  
    while(getline(infile, line))
      {
        if(line[0]=='>')
        {
        title=line;
        return true; 
        }
        else sequence +=line;
       }
 if(sequence!="") return true;
 return false;      
}
//--------------------------------------------------------------------------------------
ssize_t ReadFast(string prefixFile, string originalFile, string outputFile)
{
	ifstream preffile, origfile;
	openck(prefixFile, preffile, std::ios::in);
	openck(originalFile, origfile, std::ios::in);
	
          ofstream concatOut;
          openck(outputFile, concatOut, std::ios::out);
	

    string pref_sequence = "";
    string pref_title="",pref_title2="";
    string orig_sequence = "";
    string orig_title="";
    
    ssize_t count=0,nreads=0,pos,index;

    ReadOneSeq(preffile, pref_title2, pref_sequence);
    ReadOneSeq(origfile, orig_title, orig_sequence);

    while(ReadOneSeq(preffile, pref_title, pref_sequence))
    {
         pos = pref_title2.find("index",0)+6;
         if(pos<6)
          {
           cout << "prefix reads should have index of original reads in their name, e.g., index=2343" << endl;
           exit(1);
          }
          
         index = atoi(pref_title2.substr(pos).c_str());
      
      while(count<=index && ReadOneSeq(origfile, orig_title, orig_sequence))count++;
      
      ExtendReads(pref_sequence, orig_sequence);
      nreads++;
      
      concatOut << pref_title2 << endl;      
      concatOut << pref_sequence << endl;      
    
     if(count%100000==0){cout << ".";cout.flush();}
     pref_title2=pref_title;
     }

     concatOut.close();
     preffile.close();
     origfile.close();
     
     cout << endl;
     
     return nreads;
}      
//--------------------------------------------------------------------------------------
#define MIN3(x,y,z) ((x)<(y) ? ((x)<(z) ? (x) : (z)) : ((y)<(z) ? (y) : (z)))

// P is a macro to access 1 dimensional D as if it were 2 dimensional
#define P(i,j) ((i)*(m+1)+(j))

#define D(x,y) Data[(x)*(m+1)+y]

//--------------------------------------------------------------------------------------
ssize_t EditDistance(char A[], char B[], ssize_t match_penalty, ssize_t insert_penalty, ssize_t delete_penalty)
{
  ssize_t res;
	ssize_t m,n;
  n = strlen(A);
  m = strlen(B);

	//	int *Data; // TODO: 1. verify changes to size_t; 2. Get rid of malloc, use C++ memory storage
	//	Data = (int *)malloc(sizeof(int)*(n+1)*(m+1));
	//    if (!Data) {fprintf(stderr,"Unable to malloc memory\n"); exit(-1); }

	_INT_ *Data = new _INT_[(n+1)*(m+1)];

	if (!Data) {
		cout << "Unable to allocate memory in EditDistance" << endl;
		exit(-1);
	}


    // Initialize D
    D(0,0) = 0;		
    for (ssize_t i=1;i<=n;i++)  D(i,0) = i;
    for (ssize_t j=1;j<=m;j++)  D(0,j) = j;
    
    // Calculate the D array for the edit distance
    for (ssize_t i=1;i<=n;i++)
      for (ssize_t j=1;j<=m;j++) {
	D(i,j) = MIN3(
	  D(i-1,j-1) + (A[i-1]==B[j-1] ? 0 : match_penalty), // Match or change 
	  D(i,j-1) + insert_penalty,	// Insert 
	  D(i-1,j) + delete_penalty);	// Delete 
      }
    
    res = D(n,m);		// Store edit distance
		//    free(Data);
		delete[] Data;

  return res;
}
//---------------------------------------------------------------------------

