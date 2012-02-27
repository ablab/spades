/***************************************************************************
 * Title:          SimplifyGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include "SeqReader.h"
#include "compatibility.h"
//#include "IntegralTupleStatic.h"

using namespace std;

void PrintUsage() {
  std::cout << "usage: simplifyGraph graphFile graphOutFile -vertexSize V (20) [options]" 
						<< std::endl;
  std::cout << "options: " << std::endl;
	std::cout << "  -protectedEdges edgefile\n";
	std::cout << "                         Read in sequences from 'edgefile', and try not\n";
	std::cout << "                         to remove any edges if they contain sequences from this\n";
	std::cout << "                         file.\n" << std::endl;
  std::cout << "  -minComponentSize size Remove components less than 'size'" << std::endl
						<< "                         aggregate length. " << std::endl;
  std::cout << "  -minEdgeLength length  Remove edges that are less than 'length'" << std::endl
						<< "                         nucleotides.\n" << std::endl;
  std::cout << "  -removeSimpleBulges bulgeSize\n";
	std::cout << "                         Merge together bulges of size 'bulgeSize'\n";
  std::cout << "                         or less.\n\n";
  std::cout << "  -removeBulges bulgeSize First calculate Directed Minimum Spanning \n";
	std::cout << "                         Tree on the de Bruijn graph.  Next, add back edges \n";
	std::cout << "                         if they do not cause a cycle of length 'bulgeSize' in\n";
	std::cout << "                         the graph.\n" << std::endl;
	std::cout << "  -removeWhirls whirlSize Remove short directed cycles (src==dest) of length <\n"
						<< "                         whirlSize.\n" << std::endl;
  std::cout << "  -removeLowCoverage sigma absolute \n";
	std::cout << "                         Remove edges that have unexpectedly low coverage: less\n";
	std::cout << "                         than sigma standard deviations reads mapped to them.\n";
	std::cout << "                         Also remove edges that have less than 'cutoff', no matter\n";
	std::cout << "                         what length." << std::endl << endl;
	std::cout << "  -removeDisjointEdges   Cut edges off of vertices that do not have reads passing\n"
						<< "                         through them." <<std::endl << endl;
	std::cout << "  -skipIntervals         Do not try and read the maps of the reads to the " << endl
						<< "                         graph.  This is convenient when the graph is very " << endl
						<< "                         messy and much of it is discarded when simplifying." << endl
						<< "                         The reads must be mapped later on." << endl << endl
						<< "  -useDMST (default off) Retain all edges contained in a directed minimal spanning tree" << endl
						<< "                         when removing bulges (default is to try and remove short edges)." << endl
						<< "  -MST                   Only keep edges in the minimal spanning tree" << endl
						<< "  -intvOut file          Output intervals to 'file' (default graphOutFile.intv)" << endl
						<< "  -printSmallComponentReads readsfile componentreadsfile" << endl
						<< "                         print reads in small components to 'componentreadsfile'" << endl // TODO: what is 'readsfile' for?
						<< "  -removeLowPaths lowpathedges lowpathextend" << endl
						<< "                         Remove edges with fewer than 'lowpathedges' reads" << endl
						<< "                         A read counts if it is wholly contained in the end;" << endl
						<< "                         if the edge is wholly contained in the read;" << endl
						<< "                         or if the read extends at least 'lowpathextend' nucleotides into the read" << endl
						<< "  -suspectHPL n          Remove suspect edges with homopolymers of length n" << endl // TODO: verify and explain better
						<< "  -numReads n            " << endl; // TODO: what is it for?
		//						<< " -removeSuspectBulges n  remove suspect bulges with homopolymers of length n" << endl
		
	
  std::cout << "simplifyGraph - read in a condensed de Bruijn graph " << std::endl
	    << "  and perform bulge/whirl removal, and erosion, " << std::endl;
}

int main(int argc, char* argv[]) {

  IntervalGraph graph;
  std::string baseInName, baseOutName, edgeOutName, graphOutName, edgeFileName, pathFileName;
  std::string graphFileName, intervalFileName, bGraphOutName, intvOutName, pathOutName,
    gvzOutName;

	std::string readsFile, componentReadsFile, reportFileName;
	std::string protectedEdges;
	protectedEdges = "";
	readsFile = "";
	componentReadsFile = "";
  ssize_t minComponentSize = 0;
  ssize_t minEdgeLength = 0;
  ssize_t removeSimpleBulges = 0;
  ssize_t removeSuspectBulges = 0;
  ssize_t bulgeLength = 0;
  ssize_t suspectHPL = 0;
	ssize_t numReads = 0;
  int argi = 1;
	ssize_t removeLowPathEdges = 0;
	ssize_t lowPathEdges = 0;
	ssize_t lowPathExtend = 0;
	//UNUSED// ssize_t minSpanningReads = 2;
	ssize_t removeDisjointEdges = 0;
	double lowCoverageStddev = 0;
	ssize_t whirlLength = 0;
  if (argc < 3) {
    PrintUsage();
    exit(1);
  }
  baseInName       = argv[argi++];
  baseOutName      = argv[argi++];
  graphFileName    = baseInName  + ".bgraph";
  intervalFileName = baseInName  + ".intv";
  edgeFileName     = baseInName  + ".edge";
	pathFileName     = baseInName  + ".path";
  intvOutName      = baseOutName + ".intv";
  bGraphOutName    = baseOutName + ".bgraph";
  graphOutName     = baseOutName + ".graph";
  edgeOutName      = baseOutName + ".edge";
  gvzOutName       = baseOutName + ".dot";
	pathOutName      = baseOutName + ".path";
	reportFileName   = FormReportName(baseOutName);
	ssize_t absoluteCutoff = 2;
  int vertexSize = 20;
  ssize_t computeMST = 0;
	_INT_ skipIntervals = 0;
	ssize_t useDMST = 0;
  while(argi < argc) {
    if (strcmp(argv[argi], "-intvOut") == 0) {
      intvOutName = argv[++argi];
    }
    else if (strcmp(argv[argi], "-minComponentSize") == 0) {
      minComponentSize = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-minEdgeLength") == 0) {
      minEdgeLength = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-removeSimpleBulges") == 0) {
      removeSimpleBulges = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-vertexSize") == 0) {
      vertexSize = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-MST") == 0) {
      computeMST = 1;
    }
		else if (strcmp(argv[argi], "-skipIntervals") == 0) {
			skipIntervals = 1;
			graph.containsIntervals = 0;
		}
    else if (strcmp(argv[argi], "-removeLowCoverage") == 0) {
			if (argi < argc-2) {
				lowCoverageStddev = atof(argv[++argi]);
				absoluteCutoff    = atoi(argv[++argi]);
			}
			else {
				PrintUsage();
				exit(1);
			}
    }
    else if (strcmp(argv[argi], "-suspectHPL") == 0) {
      suspectHPL = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-removeSuspectBulges") == 0) {
      removeSuspectBulges = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-removeBulges") == 0) {
      bulgeLength = atoi(argv[++argi]);
    }
		else if (strcmp(argv[argi], "-removeWhirls") == 0) {
			whirlLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-printSmallComponentReads") == 0) {
			readsFile          = argv[++argi];
			componentReadsFile = argv[++argi];
		}
		else if (strcmp(argv[argi], "-numReads") == 0) {
			numReads = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-protectedEdges") == 0) {
			protectedEdges = argv[++argi];
		}
		else if (strcmp(argv[argi], "-removeLowPaths") == 0) {
			removeLowPathEdges = 1;
			lowPathEdges = atoi(argv[++argi]);
			lowPathExtend = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-removeDisjointEdges") == 0) {
			removeDisjointEdges = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-useDMST") == 0) {
			useDMST = 1;
		}
		else {
      PrintUsage();
      std::cout << "Bad option " << argv[argi] << std::endl;
      exit(1);
    }
    ++argi;
  }
  graph.vertexSize = vertexSize;

	std::ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);
  // read in graph components
	std::cout << CurTimeString() << ": reading interval graph" << std::endl;
  ReadIntervalGraph( baseInName, graph, vertexSize, skipIntervals, report);
	std::cout << CurTimeString() << ": done reading interval graph" << std::endl;

	report << "vertex:\t" << graph.vertexSize << std::endl;
	ssize_t e;
	for (e = 0; e < graph.edges.size(); e++ ){ 
		graph.edges[e].index = e;
	}
	if (numReads != 0)
		graph.SetMaxReadIndex(numReads);
	//  ReadSequences(edgeFileName, graph.edges, report); // TODO: check if this is redundant with ReadIntervalGraph above

  if (removeLowPathEdges) {
    std::cout << CurTimeString() << ": removing low path edges (lowPathEdges=" << lowPathEdges << ", lowPathExtend=" << lowPathExtend << ")" << std::endl;
		graph.RemoveLowPathEdges(lowPathEdges, lowPathExtend);
	}
	if (protectedEdges != "") {
    std::cout << CurTimeString() << ": protecting edges (protectedEdges=" << protectedEdges << ", vertexsize=" << graph.vertexSize+1 << ")" << std::endl;
		graph.ProtectEdges(protectedEdges, graph.vertexSize+1, report);
	}
  if (removeSuspectBulges > 0) {
    std::cout << CurTimeString() << ": marking suspect edges to remove suspect bulges (removeSuspectBulges=" << removeSuspectBulges << ")" << std::endl;
    graph.MarkSuspectEdges(removeSuspectBulges);
	}
  
  if (suspectHPL > 0) {
		std::cout << CurTimeString() << ": marking suspect edges to remove suspect homopolymers (suspectHPL=" << suspectHPL << ")" << std::endl;
    graph.MarkSuspectEdges(suspectHPL);
	}


  if (minComponentSize > 0) {
    std::cout << CurTimeString() << ": removing components smaller than " << minComponentSize << std::endl;
    graph.RemoveSmallComponents(minComponentSize, readsFile, componentReadsFile, report);
  }

  if (minEdgeLength > 0) {
    std::cout << CurTimeString() << ": eroding edges less than " << minEdgeLength << std::endl;
    graph.Erode(minEdgeLength);
  }

  if (lowCoverageStddev > 0 ) {
    std::cout << CurTimeString() << ": Removing low coverage edges " << std::endl;
    graph.RemoveLowCoverageEdges(lowCoverageStddev, absoluteCutoff);
  }
  if (removeSimpleBulges > 0) {
    std::cout << CurTimeString() << ": removing simple bulges less than " << removeSimpleBulges << std::endl;
    graph.RemoveAllSimpleBulges(removeSimpleBulges);
  }
  if (bulgeLength > 0) {
    std::cout << CurTimeString() << ": removing bulges less than " << bulgeLength << std::endl;
		//		int bl;
		//		for (bl = graph.vertexSize*2+4; bl <= bulgeLength + 4; bl+= graph.vertexSize) {
		//		std::cout << CurTimeString() << ": removing bulges of size: " << bulgeLength << std::endl;
		graph.RemoveBulges(bulgeLength, useDMST);
			//		}
  }
	if (whirlLength > 0) {
		std::cout << CurTimeString() << ": removing whirls less than " << whirlLength << std::endl;
		graph.RemoveWhirls(whirlLength);
	}


	if (removeDisjointEdges) {
		std::cout << CurTimeString() << ": Removing edges that are not spanned by at least " << removeDisjointEdges << " reads." 
							<< endl;
		ssize_t numCut = 0;
		numCut = graph.CutDisjointEdges(removeDisjointEdges);
		cout << "cut a total of " << numCut << " edges." << endl;
	}
  
  if (removeSuspectBulges > 0) {
    std::cout << CurTimeString() << ": removing suspect bulges with homopolymers of length " << removeSuspectBulges << std::endl;
		std::cout << "This function has been removed. " << std::endl;
		return 0;
    // If any edge is suspect (under the 454 read error model, I need
    // to update this for other error models), mark it with the 'marked'
    // field.  This way later I can change the way I consider edges to be good or bad.
    //    graph.RemoveSuspectBulges();
		//    graph.FindAlternatePaths();
  }


  
  if (computeMST) {
    std::cout << CurTimeString() << ": just keeping the mst " << std::endl;
    graph.RemoveAllButMST();
  }

	
  std::cout << CurTimeString() << ": clean up and sanity checks" << std::endl;
	graph.RemoveTruncatedPathIntervals();
	graph.CondenseSimplePaths();

	/* If any paths were not re-routed, get rid of them. */
	graph.DiscardGappedPaths();
	graph.CheckAllPathsContinuity(1);
	assert(graph.CheckAllPathsBalance(1));
	assert(graph.CheckBalance());
	//	graph.SetMultiplicities();
  std::cout << CurTimeString() << ": printing the graph to " << graphOutName << " and " << intvOutName << std::endl;
	//  graph.PrintIntervalGraph(bGraphOutName, intvOutName, report);
	//  PrintGraph(graph.vertices,graph.edges, graphOutName, report);
	//  PrintEdges(graph.vertices, graph.edges, edgeOutName, report);
	//  GVZPrintBGraph(graph.vertices, graph.edges, gvzOutName, report);
	//	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths, report);
	WriteIntervalGraph(baseOutName, graph, graph.vertexSize, report);
	EndReport(report);
	report.close();

	graph.Free();
  std::cout << CurTimeString() << ": Done" << std::endl;

  return 0;
}
