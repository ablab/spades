//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "version-gen.hpp"
#include "version.hpp"

namespace version {

const char *major() {
    return SPADES_VERSION_MAJOR;
}

const char *minor() {
    return SPADES_VERSION_MINOR;
}

const char *patch() {
    return SPADES_VERSION_PATCH;
}

const char *flavour() {
    return SPADES_FLAVOUR;
}

const char *package() {
    return SPADES_PACKAGE_VERSION;
}

const char *refspec() {
    return SPADES_GIT_REFSPEC;
}

const char *gitrev() {
    return SPADES_GIT_SHA1;
}

}

