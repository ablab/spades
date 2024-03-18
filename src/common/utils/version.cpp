//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
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

