#include "common.hpp"

LOGGER("p.common");

string parsed_reads;
string parsed_k_l_mers;
string parsed_k_sequence;
string parsed_l_k_mers;
string parsed_k_mers;
string distance_type;
string error_log;
string parsed_l_mers;
string graph_file;
string graph2;
string threaded_graph;
string folder;
int k = 0;
int l = 0;
int readLength = 0;
int insertLength = 0;
int minIntersect = l - 1;
int inClusterMaxShift = 1;
int useKmersVertices = 0;


int fictiveSecondReads = 0;
int needPairs = 1;
int needLmers = 1;
int needSequences = 1;
int needGraph = 1;
int useExpandDefinite = 1;
int useExtractDefinite = 1;
int needRevertedPairs = 0;
int useTraceReads = 1;
int useProcessLower = 1;

/*
 * Method adds nucleotide to the side of kMer defined by direction
 */
ll pushNucleotide(ll kMer, int length, int direction, int nucl) {
	assert(direction == LEFT || direction == RIGHT );
	if (direction == RIGHT) {
		return (ll) nucl | (kMer << (2));
	} else {
		return (ll) nucl << (2 * length) | (kMer &(((ll) 1<< (2 * length))-1));
	}
}

ll popNucleotide(ll kMer, int length, int direction) {
	if (direction == RIGHT) {
		return kMer >> 2;
	} else {
		return (((ll) 1 << (2 * length - 2)) - 1) &&  kMer;
	}
}



void initConstants(string ini_file) {
	char tmp[200];
	INFO("Trying to init constants...");

	//	string folder = string("data/");
	FILE* ini = fopen(ini_file.c_str(), "r");

	assert(fscanf(ini, "Run:\n") == 0);

	assert(fscanf(ini, "needPairs = %d\n", &needPairs) == 1);
	assert(fscanf(ini, "needLmers = %d\n", &needLmers) == 1);
	assert(fscanf(ini, "needSequences = %d\n", &needSequences) == 1);
	assert(fscanf(ini, "needGraph = %d\n", &needGraph) == 1);
	assert(fscanf(ini, "useExpandDefinite = %d\n", &useExpandDefinite) == 1);
	assert(fscanf(ini, "useTraceReads = %d\n", &useTraceReads) == 1);
	assert(fscanf(ini, "useProcessLower = %d\n", &useProcessLower) == 1);
	assert(fscanf(ini, "useExtractDefinite = %d\n", &useExtractDefinite) == 1);



	assert(fscanf(ini, "Options:\n") == 0);

	assert(fscanf(ini, "k = %d\n", &k) == 1);
	assert(fscanf(ini, "l = %d\n", &l) == 1);
	assert(fscanf(ini, "readLength = %d\n", &readLength) == 1);
	assert(fscanf(ini, "insertLength = %d\n", &insertLength) == 1);
	assert(fscanf(ini, "distance_type = %s\n" , tmp) == 1);
	distance_type = string(tmp);
	assert(fscanf(ini, "fictiveSecondReads = %d\n", &fictiveSecondReads) == 1);
	assert(fscanf(ini, "needRevertedPairs = %d\n", &needRevertedPairs) == 1);
	assert(fscanf(ini, "inClusterMaxShift = %d\n", &inClusterMaxShift) == 1);
	assert(fscanf(ini, "useKmersVertices = %d\n", &useKmersVertices) == 1);

	assert(fscanf(ini, "Filenames:\n") == 0);
	assert(fscanf(ini, "work_folder = %s\n" , tmp) == 1);
	folder = string(tmp) + '/';
	char topr[20];
	sprintf(topr, "_%d_%d", readLength, insertLength);
	string suff(topr);
	string d_desc =  suff;
	sprintf(topr, "_%d_%d",k, l);
	suff += topr;
	assert(fscanf(ini, "parsed_reads = %s\n" , tmp) == 1);
	parsed_reads = folder + string(tmp) + d_desc + ".txt";
	assert(fscanf(ini, "parsed_k_l_mers = %s\n" , tmp) == 1);
	parsed_k_l_mers = folder + string(tmp) + suff + ".txt";
	assert(fscanf(ini, "parsed_l_mers = %s\n" , tmp) == 1);
	parsed_l_mers = folder + string(tmp) + suff + ".txt";

	assert(fscanf(ini, "parsed_l_k_mers = %s\n" , tmp) == 1);
	parsed_l_k_mers = folder + string(tmp) + suff + ".txt";
	assert(fscanf(ini, "parsed_k_mers = %s\n" , tmp) == 1);
	parsed_k_mers = folder + string(tmp) + suff + ".txt";

	assert(fscanf(ini, "parsed_k_sequence = %s\n" , tmp) == 1);
	parsed_k_sequence = folder + string(tmp) + suff + ".txt";
	DEBUG(parsed_k_sequence);
	assert(fscanf(ini, "compressed_graph = %s\n" , tmp) == 1);
	graph_file = folder + string(tmp) + suff + ".dot";
	assert(fscanf(ini, "intermediate_graph = %s\n" , tmp) == 1);
	graph2 = folder + string(tmp) + suff + ".dot";
	assert(fscanf(ini, "error_log = %s\n" , tmp) == 1);
	error_log = folder + string(tmp);
	minIntersect = l - 1;

	//assert()
}



/*const string parsed_reads = "data/filtered_reads";
const string parsed_k_l_mers = string("data/klmers") + suffix  + ".txt";
const string parsed_k_sequence = string("data/vertices") + suffix  + ".txt";
const string error_log = "data/error.log";
const string parsed_l_mers = string("data/lmers")  + suffix  + ".txt";
const string graph_file = string("data/grapht")  + suffix  + ".dot";
const string graph2 = string("data/graph2") + suffix  + ".dot";
const string threaded_graph = string("data/threaded_graph") +  + ".dot";
*/
