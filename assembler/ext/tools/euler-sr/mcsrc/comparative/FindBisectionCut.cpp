/***************************************************************************
 * Title:          FindBisectionCut.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <map>
#include <sstream>

// 3rd party
#include <boost/property_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/edmunds_karp_max_flow.hpp>

// my stuff
#include "GraphDefinitions.h"
#include "utils.h"
#include "InversionBins.h"


void RemoveLowDegree(AlignmentGraph &graph,
			  ssize_t boundary,
			  std::vector<std::string> &removedVertices);

ssize_t ContractGraph(AlignmentGraph &graph,
		  ssize_t boundary,
		  std::string &removedVertex);

void RemoveMinDegree(StringVector &species,
		     InversionList &invList, 
		     ssize_t minDegree);

void RemoveVertex(AlignmentGraph &graph,
		  AlignmentVertex &vertex);

void FindMinDegree(AlignmentGraph &graph,
		   AlignmentVertex &minDegreeVertex,
		   ssize_t &minDegree,
		   ssize_t &minDegreeIndex);

void FindMinDegree(InversionList &invList, 
		   ssize_t minDegree,
		   std::vector<ssize_t> &removeList);

void BuildConnectionGraph(StringVector &species, 
			  BinMap &bins,
			  AlignmentGraph &graph);

void BuildGraph(StringVector &species, 
		BinMap &bins,
		AlignmentGraph &graph);

int main(int argc, char* argv[]) {
  
  if (argc < 2) {
    std::cout << "usage: biscut binfile minBisection [minDegree]" << std::endl;
    exit(1);
  }

  std::string binFileName, binOutName;
  int argi = 2;
  ssize_t minDegreeThreshold, minBisection;
  argi = 1;

  binFileName = argv[argi++];
  //  minBisection = atoi(argv[argi++]);
  minDegreeThreshold = 0;
  if (argi < argc) {
    minDegreeThreshold = atoi(argv[argi++]);
  }
  binOutName = "";
  if (argi < argc) {
    binOutName = argv[argi++];
  }

  StringVector species;
  std::vector<ssize_t> startPos, endPos;
  BinMap binnedInversions;
  ReadBinFile(binFileName, species, binnedInversions);
  
  BinMap::iterator binIt;
  ssize_t binNumber = 0;
  AlignmentGraph *graph;
  std::ofstream graphOut;
  std::string rejectedName = "rejectedOut";
  ssize_t rejectedGraphs = 0;
  std::stringstream sstr;
  std::vector<ssize_t> discardedChars;
  ssize_t startPosIndex = 0;
  ssize_t numMaxFlowOneVertex = 0;
  ssize_t numMaxFlowNotOne = 0;
  ssize_t numPassed = 0;
  ssize_t numFailed = 0;
  AlignmentVertexIterator vertexIt, vertexEnd;
  ssize_t numVertices;
  InversionList::iterator invIt;
  /*  for (binIt = binnedInversions.begin();
      binIt != binnedInversions.end(); ++binNumber) {*/
    graph = new AlignmentGraph;
    // Change the character/species representation into a boost 
    // graph.
    BuildConnectionGraph(species, binnedInversions, *graph);

    // Get access to all the information in the graph

    VertexNamePM vertexName = get(boost::vertex_name, *graph);
    EdgeNamePM edgeName = get(boost::edge_name, *graph);
    EdgeCapacityPM capacity = get(boost::edge_capacity, *graph);
    EdgeResidualPM residual = get(boost::edge_residual_capacity, *graph);
    EdgeReversePM  reverse  = get(boost::edge_reverse, *graph);

    std::string outName = binFileName + ".dot";
    openck(outName, graphOut);
    std::cout << "writing to " << outName << std::endl;
    write_graphviz(graphOut, *graph, make_label_writer(vertexName), make_label_writer(capacity));
    graphOut.close();
    exit(0);
    
    // Remove lowest degree vertices until almost a clique is left
    std::vector<std::string> removedVertices;
    RemoveLowDegree(*graph, 2, removedVertices);

    // Make the removed vertices unknown characters.
    std::vector<std::string>::iterator rvIt, rvEnd;
    //    std::cout << "making low degree vertices unknowns " << std::endl;
    for (rvIt = removedVertices.begin(); rvIt != removedVertices.end(); ++rvIt) {
      ssize_t refIndex, qryIndex;
      //      std::cout << "finding species " << *rvIt << std::endl;
      refIndex = FindSpecies((*binIt).second, *rvIt);
      //      std::cout << "got index for " << *rvIt << " " << refIndex <<std::endl;

      if (refIndex >= 0) {
	// Since the list of reference species is a subset of the total species, only 
	// try and remove it if it was found.
	invIt = (*binIt).second.begin() + refIndex;
      
	// Get rid of this row
	for (qryIndex = 0; qryIndex < (*binIt).second[refIndex]->inversions.size(); ++qryIndex) {
	  (*binIt).second[refIndex]->inversions[qryIndex] = 2;
	}
      }
      // Get rid of this column
      qryIndex = FindSpecies(species, *rvIt);
      for (invIt = (*binIt).second.begin(); invIt != (*binIt).second.end(); ++invIt) {
	(*invIt)->inversions[qryIndex] = 2;
      }
    }

    /*
      // This is the max flow code that's now removed. 

    // Find out what edges are left.

    boost::tie(vertexIt, vertexEnd) = boost::vertices(*graph);
    ListTraits::vertex_descriptor source, sink;
    
      std::vector<boost::default_color_type> color(boost::num_vertices(*graph));
      std::vector<ListTraits::edge_descriptor> pred(boost::num_vertices(*graph));
      //    std::vector<AlignmentEdge> pred(boost::num_vertices(*graph));
      std::cout << "operating on graph of size: " << numVertices << std::endl;
      ssize_t minDegree = -1;
      ssize_t minDegreeIndex  = -1;
      ssize_t index = 0;
      if (numVertices > 0) {
      source = *vertexIt;
      long flow;
      ssize_t vertexInd = 1;
      for ( ++vertexIt; vertexIt != vertexEnd; ++vertexIt) {
      sink = *vertexIt;
      flow = boost::edmunds_karp_max_flow(*graph, source, sink, 
      capacity, residual, reverse, &color[0], &pred[0] );
      std::cout << "flow: " << flow << " " <<  vertexName[source] << " -> " 
      << vertexName[sink] << std::endl;
      }
      // For posterity, check the minimum degree vertex
      std::cout << "min degree vertex: " << vertexName[minDegreeVertex] << " " << minDegree << std::endl;
      std::cout << "------------------------------------------------------------" << std::endl;
      if (minDegree == (ssize_t) flow) {
      numMaxFlowOneVertex++;
      }
      else {
	numMaxFlowNotOne++;
	}
	}
    if (minDegree < minBisection or numVertices <= 0) {
      std::string outName = "";
      sstr.str(outName);
      std::cout << "removing " << binNumber << std::endl;
      sstr << "rejected_" << binNumber << ".dot";
      outName = sstr.str();
      openck(outName, graphOut);
      write_graphviz(graphOut, *graph, make_label_writer(vertexName), make_label_writer(edgeName));
      graphOut.close();
      discardedChars.push_back(binNumber);
      BinMap::iterator prevBin;
      prevBin = binIt;
      ++binIt;
      binnedInversions.erase(binNumber);
      startPos.erase(startPos.begin() + startPosIndex);
      endPos.erase(endPos.begin() + startPosIndex);
      numFailed++;
    }
    else {
      ++binIt;
      ++startPosIndex;
      numPassed++;
    }
    */
    /*
      ++binIt;
      }
    */
  ssize_t b;
  std::cout << "chars that failed the minimum bisection test: " << std::endl;
  for (b = 0; b < discardedChars.size(); b++) {
    std::cout << discardedChars[b] << std::endl;
  }
  if (binOutName != "") {
    std::ofstream binOut;
    openck(binOutName, binOut);
    PrintBins(species, binnedInversions, binOut);
    binOut.close();
  }
  std::cout << "stats: " << numPassed << " passed " << numFailed << " failed " << minBisection 
	    << " singletons: " << numMaxFlowOneVertex << " other: " << numMaxFlowNotOne << std::endl;
  return 0;
}

void RemoveMinDegree(StringVector &species,
		     InversionList &invList, 
		     ssize_t minDegree) {
  std::vector<ssize_t> removeList;
  InversionList::iterator invIt, removeIt;

  FindMinDegree(invList, minDegree, removeList);

  ssize_t refSpec, qrySpec;
  ssize_t refIndex;
  refIndex = 0;
  for (invIt = invList.begin(); 
       invIt != invList.end() and removeList.size() > 0; 
       ++invIt, ++refIndex) {
    if (removeList[0] == refIndex) {
      // pop this element
      removeList.erase(removeList.begin());
      // Remove this species from the inversion list (make it unknown).
      // (remove this row)
      //      std::cout << "erasing " << (*invIt)->species << " ";
     
      for (qrySpec = 0; qrySpec < (*invIt)->inversions.size(); ++qrySpec) {
	std::cout << (*invIt)->inversions[qrySpec] << " ";
	(*invIt)->inversions[qrySpec] = 2;
      }
      std::cout << std::endl;
      // The list should be symmetric, remove the column 
      // corresponding to this species.
      for (refSpec = 0; refSpec < species.size(); ++refSpec) 
	if (species[refSpec] == (*invIt)->species)
	  break;
      //      std::cout << "erasing column: " << species[refSpec] ;
      for (removeIt =  invList.begin(); 
	   removeIt != invList.end(); ++removeIt) {
	//	std::cout << " " << (*removeIt)->inversions[refSpec];
	(*removeIt)->inversions[refSpec] = 2;
      }
      //      std::cout << std::endl;
    }
  }
}


void FindMinDegree(InversionList &invList, 
		   ssize_t minDegree,
		   std::vector<ssize_t> &removeList) {
  ssize_t refSpec, qrySpec;
  ssize_t inv;
  InversionList::iterator invIt;

  refSpec = 0;
  ssize_t degree, invVal;
  for (invIt = invList.begin(); invIt != invList.end(); ++invIt, ++refSpec) {
    degree = 0;
    for (qrySpec = 0; qrySpec < (*invIt)->inversions.size(); qrySpec++) {
      if (refSpec == qrySpec)
	continue;
      invVal = (*invIt)->inversions[qrySpec];
      //      std::cout << "invval: " << invVal << std::endl;
      if (invVal == 0 or invVal == 1)
	++degree;
    }
    //    std::cout << (*invIt)->species << " degree " << degree << " " << minDegree << std::endl;
    if (degree < minDegree)
      removeList.push_back(refSpec);
  }
}

void FindMinDegree(AlignmentGraph &graph,
		   AlignmentVertex &minDegreeVertex,
		   ssize_t &minDegree,
		   ssize_t &minDegreeIndex) {
  AlignmentVertexIterator vertexIt, vertexEnd;
  minDegree = -1;
  ssize_t index = 0;  
  
  VertexNamePM vertexName = get(boost::vertex_name, graph);
  for (boost::tie(vertexIt, vertexEnd) = boost::vertices(graph);
       vertexIt != vertexEnd;
       ++vertexIt, ++index) {
    if (minDegree < 0 or minDegree > boost::out_degree(*vertexIt, graph)) {
      minDegree = boost::out_degree(*vertexIt, graph);
      minDegreeVertex = *vertexIt;
      minDegreeIndex = index;
    }
  }
  /*
    if (minDegree >= 0) {
    std::cout << "got min degree vertex " << vertexName[minDegreeVertex] << std::endl;
  }
  */
}

void RemoveVertex(AlignmentGraph &graph,
		  AlignmentVertex &vertex) {
  boost::clear_vertex(vertex, graph);
  boost::remove_vertex(vertex, graph);
}


void RemoveLowDegree(AlignmentGraph &graph,
		     ssize_t boundary, std::vector<std::string> &removedVertices) {
  std::string removedVertex;
  while (ContractGraph(graph, boundary, removedVertex)) {
    //    std::cout << "removing species '" << removedVertex << "'" << std::endl;
    removedVertices.push_back(removedVertex);
  }
}

ssize_t ContractGraph(AlignmentGraph &graph,
		  ssize_t boundary,
		  std::string &removedVertex) {
  AlignmentVertex minDegreeVertex;
  ssize_t graphSize, minDegree, minDegreeIndex;

  graphSize = boost::num_vertices(graph);
  
  FindMinDegree(graph, minDegreeVertex, minDegree, minDegreeIndex);
  
  if (minDegree + boundary < graphSize) {
    VertexNamePM vertexName = get(boost::vertex_name, graph);
    removedVertex = vertexName[minDegreeVertex];
    RemoveVertex(graph, minDegreeVertex );
    return 1;
  }
  else {
    return 0;
  }
}
  

void BuildGraph(StringVector &species, 
		BinMap &bins,
		AlignmentGraph &graph) {
  ssize_t inv;
  InversionList::iterator invIt;
  EdgeVector edges;
  ssize_t refSpec, qrySpec;
  refSpec = 0;
  ssize_t invVal;
  // Step 1. Don't consider unconnected vertices.
  std::vector<ssize_t> linkFound;
  std::map<std::string, AlignmentVertex> vertexMap;
  VertexNamePM vertexName = get(boost::vertex_name, graph);
  EdgeNamePM edgeName = get(boost::edge_name, graph);
  EdgeCapacityPM edgeCapacity = get(boost::edge_capacity, graph);
  EdgeReversePM  reverse = get(boost::edge_reverse, graph);
  linkFound.resize(species.size());
  ssize_t s;
  for (s = 0; s < species.size(); s++) 
    linkFound[s] = 0;

  BinMap::iterator binIt;
  for (binIt = bins.begin(); binIt != bins.end(); ++binIt) {
    refSpec = 0;
    for (invIt = (*binIt).second.begin(); invIt != (*binIt).second.end(); ++invIt) {
      refSpec = FindSpecies(species, (*invIt)->species);
      for (qrySpec = 0; qrySpec < (*invIt)->inversions.size(); qrySpec++) {
	if (refSpec == qrySpec)
	  continue;
	
	invVal = (*invIt)->inversions[qrySpec];
	if (invVal == 1) {
	  linkFound[refSpec]++;
	linkFound[qrySpec]++;
	}
      }
    }
  }
  //  graph = new AlignmentGraph(species.size());
  AlignmentVertex vertex, refSpecVertex, qrySpecVertex;
  AlignmentEdge edge;
  ssize_t numCreatedVertices = 0;
  
  for (binIt = bins.begin(); binIt != bins.end(); ++binIt) {
    
    for (invIt = (*binIt).second.begin(); invIt != (*binIt).second.end(); ++invIt) {
      refSpec = FindSpecies(species, (*invIt)->species);
      if (linkFound[refSpec]) {
	if (vertexMap.find((*invIt)->species) == vertexMap.end()) {
	  vertex = boost::add_vertex(graph);
	  vertexName[vertex] = (*invIt)->species;
	  vertexMap[(*invIt)->species] = vertex;
	  //	std::cout << "adding vertex " << (*invIt)->species << std::endl;
	  ++numCreatedVertices;
	}
      }
      assert ((*invIt)->inversions.size() == species.size());
      for (qrySpec = 0; qrySpec < (*invIt)->inversions.size(); qrySpec++) {
	if (qrySpec == refSpec) continue;
	// It's possible no vertex has been added for the query vertex.  Try to find it.
	if (linkFound[qrySpec] and vertexMap.find(species[qrySpec]) == vertexMap.end()) {
	  vertex = boost::add_vertex(graph);
	  vertexName[vertex] = species[qrySpec];
	  vertexMap[species[qrySpec]] = vertex;
	  //	std::cout << "adding vertex " << species[qrySpec] << std::endl;
	  numCreatedVertices++;
	}
	invVal = (*invIt)->inversions[qrySpec];
	if (invVal == 1) {
	  // Look to see if this edge already exists.
	  ssize_t edgeInserted;
	  AlignmentVertex v1, v2;
	  AlignmentEdge   e1, e2;
	  ssize_t e1i, e2i;
	  v1 = vertexMap[(*invIt)->species];
	  v2 = vertexMap[species[qrySpec]];
	  
	  boost::graph_traits<AlignmentGraph>::out_edge_iterator edgeIt, edgeEnd;
	  edgeInserted = 0;
	  for (tie(edgeIt, edgeEnd) = boost::out_edges(v1, graph);
	       edgeIt != edgeEnd; ++edgeIt) {
	    if (target(*edgeIt, graph) == v2)
	      edgeInserted = 1;
	  }
	  
	  if (!edgeInserted) {
	    boost::tie(e1, e1i) = add_edge(v1, v2, graph);
	    boost::tie(e2, e2i) = add_edge(v2, v1, graph);
	    
	    edgeName[e1] = (*invIt)->inversions[qrySpec];
	    edgeCapacity[e1] = 1;
	    reverse[e1] = e2;
	    
	    edgeName[e2] = (*invIt)->inversions[qrySpec];
	    edgeCapacity[e2] = 1;
	    reverse[e2] = e1;
	  }
	}
      }
    }
  }
  std::cout << "created " << numCreatedVertices << std::endl;
  // Step 2. Build the edge list.
}


void BuildConnectionGraph(StringVector &species, 
			  BinMap &bins,
			  AlignmentGraph &graph) {
  ssize_t inv;
  InversionList::iterator invIt;
  EdgeVector edges;
  ssize_t refSpec, qrySpec;
  refSpec = 0;
  ssize_t invVal;
  // Step 1. Don't consider unconnected vertices.
  std::vector<ssize_t> linkFound;
  std::map<std::string, AlignmentVertex> vertexMap;
  VertexNamePM vertexName = get(boost::vertex_name, graph);
  EdgeNamePM edgeName = get(boost::edge_name, graph);
  EdgeCapacityPM edgeCapacity = get(boost::edge_capacity, graph);
  EdgeReversePM  reverse = get(boost::edge_reverse, graph);
  linkFound.resize(species.size());
  ssize_t s;
  for (s = 0; s < species.size(); s++) 
    linkFound[s] = 0;
  ssize_t refVal, qryVal;
  refSpec = 0;
  AlignmentVertex vertex, refSpecVertex, qrySpecVertex;
  AlignmentEdge edge;
  BinMap::iterator binIt;
  for (binIt = bins.begin(); binIt != bins.end(); ++binIt) {
    for (invIt = (*binIt).second.begin(); invIt != (*binIt).second.end(); ++invIt) {
      for (refSpec = 0; refSpec < (*invIt)->inversions.size() - 1; refSpec++) {
	for (qrySpec = refSpec + 1; qrySpec < (*invIt)->inversions.size(); qrySpec++) {
	  refVal = (*invIt)->inversions[refSpec];
	  qryVal = (*invIt)->inversions[qrySpec];
	  
	  if (refVal == 1 and qryVal == 1) {
	    linkFound[refSpec]++;
	    linkFound[qrySpec]++;
	    
	    if (vertexMap.find(species[refSpec]) == vertexMap.end()) {
	      vertex = boost::add_vertex(graph);
	    vertexName[vertex] = species[refSpec];
	    vertexMap[species[refSpec]] = vertex;
	    //	std::cout << "adding vertex " << (*invIt)->species << std::endl;
	    }
	    if (vertexMap.find(species[qrySpec]) == vertexMap.end()) {
	      vertex = boost::add_vertex(graph);
	    vertexName[vertex] = species[qrySpec];
	    vertexMap[species[qrySpec]] = vertex;
	    //	std::cout << "adding vertex " << species[qrySpec] << std::endl;
	    }
	    
	    // Look to see if this edge already exists.
	    ssize_t edgeInserted;
	    AlignmentVertex v1, v2;
	    AlignmentEdge   e1, e2;
	    ssize_t e1i, e2i;
	  v1 = vertexMap[species[refSpec]];
	  v2 = vertexMap[species[qrySpec]];
	  
	  boost::graph_traits<AlignmentGraph>::out_edge_iterator edgeIt, edgeEnd;
	  edgeInserted = 0;
	  for (tie(edgeIt, edgeEnd) = boost::out_edges(v1, graph);
	       edgeIt != edgeEnd and ! edgeInserted; ++edgeIt) {
	    if (target(*edgeIt, graph) == v2) {
	      //	      std::cout << "already added an edge " << std::endl;
	      edgeInserted = 1;
	      edgeCapacity[*edgeIt]++;
	      //	      std::cout << "cpacity: " << edgeCapacity[*edgeIt] << std::endl;
	    }
	  }
	  
	  if (!edgeInserted) {
	    boost::tie(e1, e1i) = add_edge(v1, v2, graph);
	    //	    boost::tie(e2, e2i) = add_edge(v2, v1, graph);
	    
	    edgeName[e1] = (*invIt)->inversions[qrySpec];
	    edgeCapacity[e1] = 1;
	    //	    reverse[e1] = e2;
	    /*	    
	      edgeName[e2] = (*invIt)->inversions[qrySpec];
	    edgeCapacity[e2] = 1;
	    reverse[e2] = e1;
	    */
	  }
	  }
	}
      }
    }
  }
  // Step 2. Build the edge list.
}
