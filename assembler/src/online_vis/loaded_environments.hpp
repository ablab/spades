//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "environment.hpp"

namespace online_visualization {

    template <class Env>
    class LoadedEnvironments : public map<string, shared_ptr<Env> > {};
//    typedef map<string, shared_ptr<Env> > LoadedEnvironments;
}
