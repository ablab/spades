#include "common.hpp"

LOGGER("p.common");

string parsed_reads;
string parsed_k_l_mers;
string parsed_k_sequence;
string error_log;
string parsed_l_mers;
string graph_file;
string graph2;
string threaded_graph;

int k = 0;
int l = 0;
int readLength = 0;
int insertLength = 0;
int minIntersect;


int needPairs = 1;
int needLmers = 1;
int needSequences = 1;
int needGraph = 1;
int useExpandDefinite = 1;
int useTraceReads = 1;
int useProcessLower = 1;
void initConstants(string ini_file) {
	char tmp[200];
	INFO("Trying to init constants...");

	//	string folder = string("data/");
	FILE* ini = fopen(ini_file.c_str(), "r");
	assert(fscanf(ini, "k = %d\n", &k) == 1);
	assert(fscanf(ini, "l = %d\n", &l) == 1);
	assert(fscanf(ini, "readLength = %d\n", &readLength) == 1);
	assert(fscanf(ini, "insertLength = %d\n", &insertLength) == 1);
//	assert(fscanf(ini, "maxSeqLength = %d\n", &maxSeqLength) == 1);
	assert(fscanf(ini, "Filenames:\n") == 0);
	assert(fscanf(ini, "work_folder = %s\n" , tmp) == 1);
	string folder = string(tmp) + '/';
	assert(fscanf(ini, "distance_type = %s\n" , tmp) == 1);
	string suff("_");
	suff += readLength;
	suff += "_" ;
	suff += insertLength;
	suff += "_";
	suff += string(tmp);
	char topr[20];
	sprintf(topr, "_%d_%d",k, l);
	suff += topr;
	string d_desc = "_" + string(tmp);
	assert(fscanf(ini, "parsed_reads = %s\n" , tmp) == 1);
	parsed_reads = folder + string(tmp) + d_desc + ".txt";
	assert(fscanf(ini, "parsed_k_l_mers = %s\n" , tmp) == 1);
	parsed_k_l_mers = folder + string(tmp) + suff + ".txt";
	assert(fscanf(ini, "parsed_l_mers = %s\n" , tmp) == 1);
	parsed_l_mers = folder + string(tmp) + suff + ".txt";
	assert(fscanf(ini, "parsed_k_sequence = %s\n" , tmp) == 1);
	parsed_k_sequence = folder + string(tmp) + suff + ".txt";
	assert(fscanf(ini, "compressed_graph = %s\n" , tmp) == 1);
	graph_file = folder + string(tmp) + suff + ".dot";
	assert(fscanf(ini, "intermediate_graph = %s\n" , tmp) == 1);
	graph2 = folder + string(tmp) + suff + ".dot";
	assert(fscanf(ini, "error_log = %s\n" , tmp) == 1);
	error_log = folder + string(tmp);
	assert(fscanf(ini, "Run:\n") == 0);

	assert(fscanf(ini, "needPairs = %d\n", &needPairs) == 1);
	assert(fscanf(ini, "needLmers = %d\n", &needLmers) == 1);
	assert(fscanf(ini, "needSequences = %d\n", &needSequences) == 1);
	assert(fscanf(ini, "needGraph = %d\n", &needGraph) == 1);
	assert(fscanf(ini, "useExpandDefinite = %d\n", &useExpandDefinite) == 1);
	assert(fscanf(ini, "useTraceReads = %d\n", &useTraceReads) == 1);
	assert(fscanf(ini, "useProcessLower = %d\n", &useProcessLower) == 1);
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
