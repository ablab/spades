//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef GAP_CLOSER_HPP_
#define GAP_CLOSER_HPP_

#include "stage.hpp"

namespace debruijn_graph {

class GapClosing : public spades::AssemblyStage {
  public:
    GapClosing(const char* id)
        : AssemblyStage("Gap Closer", id) {}

    void run(conj_graph_pack &gp);
};

}



#endif /* GAP_CLOSER_HPP_ */
