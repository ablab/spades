/***************************************************************************
 * Title:          GraphDefinitions.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef GRAPH_DEFINITIONS_H_
#define GRAPH_DEFINITIONS_H_

#include <string>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

typedef std::pair<ssize_t, ssize_t> Edge;
typedef std::vector<Edge> EdgeVector;

// This is what caused me some problems.  The first two arguments are the 
// storage types for vertices and edges.  They must be the same as those used in the 
// graph definition.  I had mapS and mapS for this, and it wouldn't compile.  It seems
// obvious in retrospect.

typedef boost::adjacency_list_traits<boost::mapS, boost::vecS, boost::directedS> ListTraits;

typedef boost::property<boost::vertex_name_t, std::string> VertexProperty;

typedef boost::property<boost::edge_name_t, ssize_t,
         boost::property<boost::edge_capacity_t, ssize_t,
          boost::property<boost::edge_residual_capacity_t, ssize_t,
	   boost::property<boost::edge_reverse_t, ListTraits::edge_descriptor> 
          > 
         > 
        > EdgeProperty;

typedef boost::adjacency_list<boost::mapS, boost::vecS,
			      boost::directedS,
			      VertexProperty, EdgeProperty> AlignmentGraph;
// Convenient to make a typedef for.
typedef boost::graph_traits<AlignmentGraph>::vertex_descriptor AlignmentVertex;
typedef boost::graph_traits<AlignmentGraph>::edge_descriptor   AlignmentEdge;
typedef boost::graph_traits<AlignmentGraph>::vertex_iterator   AlignmentVertexIterator;
typedef boost::property_map<AlignmentGraph, boost::vertex_name_t>::type VertexNamePM;
typedef boost::property_map<AlignmentGraph, boost::edge_name_t>::type EdgeNamePM;
typedef boost::property_map<AlignmentGraph, boost::edge_capacity_t>::type EdgeCapacityPM;
typedef boost::property_map<AlignmentGraph, boost::edge_residual_capacity_t>::type EdgeResidualPM;
typedef boost::property_map<AlignmentGraph, boost::edge_reverse_t>::type EdgeReversePM;

#endif
