/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "run.h"

static void printUsage()
{
	puts("Usage:");
	puts("./velvetg directory [options]");
	puts("");
	puts("\tdirectory\t\t\t: working directory name");
	puts("");
	puts("Standard options:");
	puts("\t-cov_cutoff <floating-point|auto>\t: removal of low coverage nodes AFTER tour bus or allow the system to infer it");
	puts("\t\t(default: no removal)");
	puts("\t-ins_length <integer>\t\t: expected distance between two paired end reads (default: no read pairing)");
	puts("\t-read_trkg <yes|no>\t\t: tracking of short read positions in assembly (default: no tracking)");
	puts("\t-min_contig_lgth <integer>\t: minimum contig length exported to contigs.fa file (default: hash length * 2)");
	puts("\t-amos_file <yes|no>\t\t: export assembly to AMOS file (default: no export)");
	puts("\t-exp_cov <floating point|auto>\t: expected coverage of unique regions or allow the system to infer it");
	puts("\t\t(default: no long or paired-end read resolution)");
	puts("\t-long_cov_cutoff <floating-point>: removal of nodes with low long-read coverage AFTER tour bus");
	puts("\t\t(default: no removal)");
	puts("");
	puts("Advanced options:");
	puts("\t-ins_length2 <integer>\t\t: expected distance between two paired-end reads in the second short-read dataset (default: no read pairing)");
	puts("\t-ins_length_long <integer>\t: expected distance between two long paired-end reads (default: no read pairing)");
	puts("\t-ins_length*_sd <integer>\t: est. standard deviation of respective dataset (default: 10% of corresponding length)");
	puts("\t\t[replace '*' by nothing, '2' or '_long' as necessary]");
	puts("\t-scaffolding <yes|no>\t\t: scaffolding of contigs used paired end information (default: on)");
	puts("\t-max_branch_length <integer>\t: maximum length in base pair of bubble (default: 100)");
	puts("\t-max_divergence <floating-point>: maximum divergence rate between two branches in a bubble (default: 0.2)");
	puts("\t-max_gap_count <integer>\t: maximum number of gaps allowed in the alignment of the two branches of a bubble (default: 3)");
	puts("\t-min_pair_count <integer>\t: minimum number of paired end connections to justify the scaffolding of two long contigs (default: 5)");
	puts("\t-max_coverage <floating point>\t: removal of high coverage nodes AFTER tour bus (default: no removal)");
	puts("\t-coverage_mask <int>\t: minimum coverage required for confident regions of contigs (default: 1)");
	puts("\t-long_mult_cutoff <int>\t\t: minimum number of long reads required to merge contigs (default: 2)");
	puts("\t-unused_reads <yes|no>\t\t: export unused reads in UnusedReads.fa file (default: no)");
	puts("\t-alignments <yes|no>\t\t: export a summary of contig alignment to the reference sequences (default: no)");
	puts("\t-exportFiltered <yes|no>\t: export the long nodes which were eliminated by the coverage filters (default: no)");
	puts("\t-clean <yes|no>\t\t\t: remove all the intermediary files which are useless for recalculation (default : no)");
	puts("\t-very_clean <yes|no>\t\t: remove all the intermediary files (no recalculation possible) (default: no)");
	puts("\t-paired_exp_fraction <double>\t: remove all the paired end connections which less than the specified fraction of the expected count (default: 0.1)");
	puts("\t-shortMatePaired* <yes|no>\t: for mate-pair libraries, indicate that the library might be contaminated with paired-end reads (default no)");
	puts("");
	puts("Output:");
	puts("\tdirectory/contigs.fa\t\t: fasta file of contigs longer than twice hash length");
	puts("\tdirectory/stats.txt\t\t: stats file (tab-spaced) useful for determining appropriate coverage cutoff");
	puts("\tdirectory/LastGraph\t\t: special formatted file with all the information on the final graph");
	puts("\tdirectory/velvet_asm.afg\t: (if requested) AMOS compatible assembly file");
}

int main(int argc, char **argv)
{
	ReadSet *sequences = NULL;
	RoadMapArray *rdmaps;
	PreGraph *preGraph;
	Graph *graph;
	char *directory, *graphFilename, *preGraphFilename, *seqFilename,
	    *roadmapFilename, *lowCovContigsFilename, *highCovContigsFilename;
	double coverageCutoff = -1;
	double longCoverageCutoff = -1;
	double maxCoverageCutoff = -1;
	double expectedCoverage = -1;
	Coordinate minContigLength = -1;
	Coordinate minContigKmerLength;
	boolean *dubious = NULL;
	Coordinate insertLength[CATEGORIES];
	Coordinate insertLengthLong = -1;
	Coordinate std_dev[CATEGORIES];
	Coordinate std_dev_long = -1;
	short int accelerationBits = 24;
	boolean readTracking = false;
	boolean exportAssembly = false;
	boolean unusedReads = false;
	boolean estimateCoverage = false;
	boolean estimateCutoff = false;
	boolean exportAlignments = false;
	FILE *file;
	int arg_index, arg_int;
	double arg_double;
	char *arg;
	ShortLength *sequenceLengths = NULL;
	Category cat;
	boolean scaffolding = true;
	int pebbleRounds = 1;
	long long longlong_var;
	short int short_var;
	boolean exportFilteredNodes = false;
	int clean = 0;
	boolean shadows[CATEGORIES];
	int coverageMask = 1;

	setProgramName("velvetg");

	for (cat = 0; cat < CATEGORIES; cat++) {
		insertLength[cat] = -1;
		std_dev[cat] = -1;
		shadows[cat] = false;
	}

	// Error message
	if (argc == 1) {
		puts("velvetg - de Bruijn graph construction, error removal and repeat resolution");
		printf("Version %i.%i.%2.2i\n", VERSION_NUMBER,
		       RELEASE_NUMBER, UPDATE_NUMBER);
		puts("Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)");
		puts("This is free software; see the source for copying conditions.  There is NO");
		puts("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.");
		puts("Compilation settings:");
		printf("CATEGORIES = %i\n", CATEGORIES);
		printf("MAXKMERLENGTH = %i\n", MAXKMERLENGTH);
#ifdef _OPENMP
		puts("OPENMP");
#endif
#ifdef LONGSEQUENCES
		puts("LONGSEQUENCES");
#endif
#ifdef BIGASSEMBLY
		puts("BIGASSEMBLY");
#endif
#ifdef COLOR
		puts("COLOR");
#endif
#ifdef DEBUG
		puts("DEBUG");
#endif
		puts("");
		printUsage();
		return 1;
	}

	if (strcmp(argv[1], "--help") == 0) {
		printUsage();
		return 0;
	}

	// Memory allocation 
	directory = argv[1];
	graphFilename = mallocOrExit(strlen(directory) + 100, char);
	preGraphFilename =
	    mallocOrExit(strlen(directory) + 100, char);
	roadmapFilename = mallocOrExit(strlen(directory) + 100, char);
	seqFilename = mallocOrExit(strlen(directory) + 100, char);
	lowCovContigsFilename = mallocOrExit(strlen(directory) + 100, char);
	highCovContigsFilename = mallocOrExit(strlen(directory) + 100, char);
	// Argument parsing
	for (arg_index = 2; arg_index < argc; arg_index++) {
		arg = argv[arg_index++];
		if (arg_index >= argc) {
			velvetLog("Unusual number of arguments!\n");
			printUsage();
#ifdef DEBUG 
			abort();
#endif 
			exit(1);
		}

		if (strcmp(arg, "-cov_cutoff") == 0) {
			if (strcmp(argv[arg_index], "auto") == 0) {
				estimateCutoff = true;
			} else {
				sscanf(argv[arg_index], "%lf", &coverageCutoff);
			}
		} else if (strcmp(arg, "-long_cov_cutoff") == 0) {
			sscanf(argv[arg_index], "%lf", &longCoverageCutoff);
		} else if (strcmp(arg, "-exp_cov") == 0) {
			if (strcmp(argv[arg_index], "auto") == 0) {
				estimateCoverage = true;
				readTracking = true;
			} else {
				sscanf(argv[arg_index], "%lf", &expectedCoverage);
				if (expectedCoverage > 0)
					readTracking = true;
			}
		} else if (strcmp(arg, "-ins_length") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			insertLength[0] = (Coordinate) longlong_var;
			if (insertLength[0] < 0) {
				velvetLog("Invalid insert length: %lli\n",
				       (long long) insertLength[0]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
		} else if (strcmp(arg, "-ins_length_sd") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			std_dev[0] = (Coordinate) longlong_var;
			if (std_dev[0] < 0) {
				velvetLog("Invalid std deviation: %lli\n",
				       (long long) std_dev[0]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
		} else if (strcmp(arg, "-ins_length_long") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			insertLengthLong = (Coordinate) longlong_var;
		} else if (strcmp(arg, "-ins_length_long_sd") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			std_dev_long = (Coordinate) longlong_var;
		} else if (strncmp(arg, "-ins_length", 11) == 0
			   && strchr(arg, 'd') == NULL) {
			sscanf(arg, "-ins_length%hi", &short_var);
			cat = (Category) short_var;
			if (cat < 1 || cat > CATEGORIES) {
				velvetLog("Unknown option: %s\n", arg);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
			sscanf(argv[arg_index], "%lli", &longlong_var);
			insertLength[cat - 1] = (Coordinate) longlong_var;
			if (insertLength[cat - 1] < 0) {
				velvetLog("Invalid insert length: %lli\n",
				       (long long) insertLength[cat - 1]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
		} else if (strncmp(arg, "-ins_length", 11) == 0) {
			sscanf(arg, "-ins_length%hi_sd", &short_var);
			cat = (Category) short_var;
			if (cat < 1 || cat > CATEGORIES) {
				velvetLog("Unknown option: %s\n", arg);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
			sscanf(argv[arg_index], "%lli", &longlong_var);
			std_dev[cat - 1] = (Coordinate) longlong_var;
			if (std_dev[cat - 1] < 0) {
				velvetLog("Invalid std deviation: %lli\n",
				       (long long) std_dev[cat - 1]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
		} else if (strcmp(arg, "-read_trkg") == 0) {
			readTracking =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-scaffolding") == 0) {
			scaffolding =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-exportFiltered") == 0) {
			exportFilteredNodes =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-amos_file") == 0) {
			exportAssembly =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-alignments") == 0) {
			exportAlignments =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-min_contig_lgth") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			minContigLength = (Coordinate) longlong_var;
		} else if (strcmp(arg, "-coverage_mask") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			coverageMask = (IDnum) longlong_var;
		} else if (strcmp(arg, "-accel_bits") == 0) {
			sscanf(argv[arg_index], "%hi", &accelerationBits);
			if (accelerationBits < 0) {
				velvetLog
				    ("Illegal acceleration parameter: %s\n",
				     argv[arg_index]);
				printUsage();
				return -1;
			}
		} else if (strcmp(arg, "-max_branch_length") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setMaxReadLength(arg_int);
			setLocalMaxReadLength(arg_int);
		} else if (strcmp(arg, "-max_divergence") == 0) {
			sscanf(argv[arg_index], "%lf", &arg_double);
			setMaxDivergence(arg_double);
			setLocalMaxDivergence(arg_double);
		} else if (strcmp(arg, "-max_gap_count") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setMaxGaps(arg_int);
			setLocalMaxGaps(arg_int);
		} else if (strcmp(arg, "-min_pair_count") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setUnreliableConnectionCutoff(arg_int);
		} else if (strcmp(arg, "-max_coverage") == 0) {
			sscanf(argv[arg_index], "%lf", &maxCoverageCutoff);
		} else if (strcmp(arg, "-long_mult_cutoff") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setMultiplicityCutoff(arg_int);
		} else if (strcmp(arg, "-paired_exp_fraction") == 0) {
			sscanf(argv[arg_index], "%lf", &arg_double);
			setPairedExpFraction(arg_double);
		} else if (strcmp(arg, "-clean") == 0) {
			if (strcmp(argv[arg_index], "yes") == 0)
				clean = 1;
		} else if (strcmp(arg, "-very_clean") == 0) {
			if (strcmp(argv[arg_index], "yes") == 0)
				clean = 2;
		} else if (strcmp(arg, "-unused_reads") == 0) {
			unusedReads =
			    (strcmp(argv[arg_index], "yes") == 0);
			if (unusedReads)
				readTracking = true;
		} else if (strcmp(arg, "-shortMatePaired") == 0) {
			shadows[0] = (strcmp(argv[arg_index], "yes") == 0);
		} else if (strncmp(arg, "-shortMatePaired", 16) == 0) {
			sscanf(arg, "-shortMatePaired%hi", &short_var);
			cat = (Category) short_var;
			if (cat < 1 || cat > CATEGORIES) {
				velvetLog("Unknown option: %s\n", arg);
#ifdef DEBUG
				abort();
#endif
				exit(1);
			}
			shadows[cat - 1] = (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "--help") == 0) {
			printUsage();
			return 0;
		} else {
			velvetLog("Unknown option: %s;\n", arg);
			printUsage();
			return 1;
		}
	}

	// Bookkeeping
	logInstructions(argc, argv, directory);

	strcpy(seqFilename, directory);
	strcat(seqFilename, "/Sequences");

	strcpy(roadmapFilename, directory);
	strcat(roadmapFilename, "/Roadmaps");

	strcpy(preGraphFilename, directory);
	strcat(preGraphFilename, "/PreGraph");

	if (!readTracking) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Graph");
	} else {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Graph2");
	}

	strcpy(lowCovContigsFilename, directory);
	strcat(lowCovContigsFilename, "/lowCoverageContigs.fa");

	strcpy(highCovContigsFilename, directory);
	strcat(highCovContigsFilename, "/highCoverageContigs.fa");

	// Graph uploading or creation
	if ((file = fopen(graphFilename, "r")) != NULL) {
		fclose(file);
		graph = importGraph(graphFilename);
	} else if ((file = fopen(preGraphFilename, "r")) != NULL) {
		fclose(file);
		sequences = importReadSet(seqFilename);
		convertSequences(sequences);
		graph =
		    importPreGraph(preGraphFilename, sequences,
				   roadmapFilename, readTracking, accelerationBits);
		sequenceLengths =
		    getSequenceLengths(sequences, getWordLength(graph));
		correctGraph(graph, sequenceLengths, sequences->categories);
		exportGraph(graphFilename, graph, sequences->tSequences);
	} else if ((file = fopen(roadmapFilename, "r")) != NULL) {
		fclose(file);

		rdmaps = importRoadMapArray(roadmapFilename);
		preGraph = newPreGraph_pg(rdmaps, seqFilename);
		concatenatePreGraph_pg(preGraph);
		clipTips_pg(preGraph);
		exportPreGraph_pg(preGraphFilename, preGraph);
		destroyPreGraph_pg(preGraph);

		sequences = importReadSet(seqFilename);
		convertSequences(sequences);
		graph =
		    importPreGraph(preGraphFilename, sequences,
				   roadmapFilename, readTracking, accelerationBits);
		sequenceLengths =
		    getSequenceLengths(sequences, getWordLength(graph));
		correctGraph(graph, sequenceLengths, sequences->categories);
		exportGraph(graphFilename, graph, sequences->tSequences);
	} else {
		velvetLog("No Roadmap file to build upon! Please run velveth (see manual)\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}

	// Set insert lengths and their standard deviations
	for (cat = 0; cat < CATEGORIES; cat++) {
		if (insertLength[cat] > -1 && std_dev[cat] < 0)
			std_dev[cat] = insertLength[cat] / 10;
		setInsertLengths(graph, cat,
				 insertLength[cat], std_dev[cat]);
	}

	if (insertLengthLong > -1 && std_dev_long < 0)
		std_dev_long = insertLengthLong / 10;
	setInsertLengths(graph, CATEGORIES,
			 insertLengthLong, std_dev_long);

	// Coverage cutoff
	if (expectedCoverage < 0 && estimateCoverage == true) {
		expectedCoverage = estimated_cov(graph, directory);
		if (coverageCutoff < 0) {
			coverageCutoff = expectedCoverage / 2;
			estimateCutoff = true;
		}
	} else { 
		estimateCoverage = false;
		if (coverageCutoff < 0 && estimateCutoff) 
			coverageCutoff = estimated_cov(graph, directory) / 2;
		else 
			estimateCutoff = false;
	}

	if (coverageCutoff < 0) {
		velvetLog("WARNING: NO COVERAGE CUTOFF PROVIDED\n");
		velvetLog("Velvet will probably leave behind many detectable errors\n");
		velvetLog("See manual for instructions on how to set the coverage cutoff parameter\n");
	}

	if (sequences == NULL) {
		sequences = importReadSet(seqFilename);
		convertSequences(sequences);
	}

	if (minContigLength < 2 * getWordLength(graph))
		minContigKmerLength = getWordLength(graph);
	else
		minContigKmerLength = minContigLength - getWordLength(graph) + 1;		

	dubious =
	    removeLowCoverageNodesAndDenounceDubiousReads(graph,
							  coverageCutoff,
							  sequences,
							  exportFilteredNodes,
							  minContigKmerLength,
							  lowCovContigsFilename);

	removeLowLongCoverageNodesAndDenounceDubiousReads(graph,
							  longCoverageCutoff,
							  sequences,
							  dubious,
							  exportFilteredNodes,
							  minContigKmerLength,
							  lowCovContigsFilename);

	removeHighCoverageNodes(graph, maxCoverageCutoff, exportFilteredNodes, minContigKmerLength, highCovContigsFilename);
	clipTipsHard(graph);

	if (sequences->readCount > 0 && sequences->categories[0] == REFERENCE)
		removeLowArcs(graph, coverageCutoff);

	if (expectedCoverage > 0) {

		// Mixed length sequencing
		readCoherentGraph(graph, isUniqueSolexa, expectedCoverage,
				  sequences);

		// Paired end resolution
		createReadPairingArray(sequences);
		pebbleRounds += pairedCategories(sequences);
		detachDubiousReads(sequences, dubious);
		activateGapMarkers(graph);
		for ( ;pebbleRounds > 0; pebbleRounds--)
			exploitShortReadPairs(graph, sequences, dubious, shadows, scaffolding);
	} else {
		velvetLog("WARNING: NO EXPECTED COVERAGE PROVIDED\n");
		velvetLog("Velvet will be unable to resolve any repeats\n");
		velvetLog("See manual for instructions on how to set the expected coverage parameter\n");
	}

	if (dubious)
		free(dubious);

	concatenateGraph(graph);

	removeLowCoverageReferenceNodes(graph, coverageCutoff, longCoverageCutoff, sequences);

	strcpy(graphFilename, directory);
	strcat(graphFilename, "/contigs.fa");
	sequenceLengths = getSequenceLengths(sequences, getWordLength(graph));
	exportLongNodeSequences(graphFilename, graph, minContigKmerLength, sequences, sequenceLengths, coverageMask); 

	if (exportAlignments) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/contig-alignments.psa");
		exportLongNodeMappings(graphFilename, graph, sequences,
					     minContigKmerLength, seqFilename);
	}

	strcpy(graphFilename, directory);
	strcat(graphFilename, "/stats.txt");
	displayGeneralStatistics(graph, graphFilename, sequences);

	if (clean == 0) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/LastGraph");
		exportGraph(graphFilename, graph, sequences->tSequences);
	}

	if (exportAssembly) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/velvet_asm.afg");
		exportAMOSContigs(graphFilename, graph, minContigKmerLength, sequences);
	}

	if (unusedReads)
		exportUnusedReads(graph, sequences, minContigKmerLength, directory);

	if (estimateCoverage) 
		velvetLog("Estimated Coverage = %f\n", expectedCoverage);
	if (estimateCutoff) 
		velvetLog("Estimated Coverage cutoff = %f\n", coverageCutoff);

	logFinalStats(graph, minContigKmerLength, directory);

	if (clean > 0) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Roadmaps");
		remove(graphFilename);	

		strcpy(graphFilename, directory);
		strcat(graphFilename, "/LastGraph");
		remove(graphFilename);	
	} 

	if (clean > 1) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Sequences");
		remove(graphFilename);	

		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Graph2.txt");
		remove(graphFilename);	

		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Graph.txt");
		remove(graphFilename);	
	}

	free(sequenceLengths);
	destroyGraph(graph);
	free(graphFilename);
	free(preGraphFilename);
	free(seqFilename);
	free(roadmapFilename);
	free(lowCovContigsFilename);
	free(highCovContigsFilename);
	destroyReadSet(sequences);
	return 0;
}
