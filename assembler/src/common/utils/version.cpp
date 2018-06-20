//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "version-gen.hpp"
#include "version.hpp"

namespace version {

const char *refspec() {
    return SPADES_GIT_REFSPEC;
}

const char *gitrev() {
    return SPADES_GIT_SHA1;
}

}

