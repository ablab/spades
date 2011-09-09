/*
 * labeler.hpp
 *
 *  Created on: 8 Sep 2011
 *      Author: valery
 */

#pragma once

#include "omni/total_labeler.hpp"
#include "omni/graph_labeler.hpp"

namespace debruijn_graph
{

class ConjugateDeBruijnGraph;

typedef omnigraph::TotalLabelerGraphStruct  <debruijn_graph::ConjugateDeBruijnGraph>	total_labeler_graph_struct;
typedef omnigraph::TotalLabeler             <debruijn_graph::ConjugateDeBruijnGraph>    total_labeler;
typedef omnigraph::GraphLabeler             <debruijn_graph::ConjugateDeBruijnGraph>    graph_labeler;

} // debruijn_graph
