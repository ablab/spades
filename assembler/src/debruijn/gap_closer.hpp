//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef GAP_CLOSER_HPP_
#define GAP_CLOSER_HPP_

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class GapClosing : public spades::AssemblyStage {
  public:
    GapClosing(const char* id)
        : AssemblyStage("Gap Closer", id) {}

    void run(conj_graph_pack &gp, const char*);
};

}



#endif /* GAP_CLOSER_HPP_ */
