//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "stage.hpp"

namespace debruijn_graph {

template<class Read>
void construct_graph(io::ReadStreamList<Read>& streams,
                     conj_graph_pack& gp, io::SingleStreamPtr contigs_stream = io::SingleStreamPtr());

class Construction : public spades::AssemblyStage {
  public:
    Construction()
        : AssemblyStage("Construction", "construction") {}

    void run(conj_graph_pack &gp, const char*);
};

}

