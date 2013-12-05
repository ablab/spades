//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * labeler.hpp
 *
 *  Created on: 8 Sep 2011
 *      Author: valery
 */

#pragma once

#include "omni/total_labeler.hpp"
#include "omni/visualization/graph_labeler.hpp"
#include "debruijn_graph.hpp"

namespace debruijn_graph
{
typedef omnigraph::TotalLabeler             <ConjugateDeBruijnGraph>    total_labeler;
typedef omnigraph::GraphLabeler             <ConjugateDeBruijnGraph>    graph_labeler;

} // debruijn_graph
